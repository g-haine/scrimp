# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Prototype: one separate binary writer process for rank-owned PETSc states.

No import of scrimp, getfem, petsc4py or mpi4py. The writer is a fresh
Python interpreter, never a forked MPI worker. Producers use short local
NumPy copies and put_nowait() into a bounded queue; only their local
background thread performs blocking socket I/O. A full queue is a
FAIL-FAST error, never an unbounded-memory buffer or silent data loss.

Linux / Unix-domain sockets, shared local filesystem. The single
writer uses one I/O thread per attached rank. An atomic manifest is
published only if ALL ranks supplied matching sample times/steps and
their owned intervals cover the entire global vector. No experiment
should treat unfinished .partial files as a valid trajectory.

Storage protocol for SCRIMP's opt-in trajectory backend. In production,
the I/O interpreter must be started separately, before MPI ranks.
"""
import argparse
import json
import os
from pathlib import Path
import queue
import socket
import socketserver
import struct
import subprocess
import sys
import threading
import time

import numpy as np


MAGIC = b"SCSTBIN1"
HEADER = struct.Struct("<8sqqqqqq")
FRAME = struct.Struct("<Bqd")
DATA, END = 1, 2
DTYPES = {1: np.dtype("<f8"), 2: np.dtype("<c16")}
CODE = {np.dtype("<f8"): 1, np.dtype("<c16"): 2}


def _recv_exact(sock, length):
    value = bytearray(length)
    view = memoryview(value)
    pos = 0
    while pos < length:
        received = sock.recv_into(view[pos:])
        if received == 0:
            raise EOFError("Writer socket closed before completing frame")
        pos += received
    return value


class AsyncOutputBackpressure(RuntimeError):
    """The bounded sender queue is full; do not discard a saved state."""


class RankSender:
    """One rank's bounded producer and private background transport thread.

    Construct and close collectively at safe setup/teardown points OUTSIDE
    the PETSc TS callback; enqueue() itself calls no MPI methods.
    """

    def __init__(self, socket_path, *, rank, ranks, global_size, lo, hi,
                 dtype=np.float64, max_queued=4, connect_timeout=10.):
        self.dtype = np.dtype(dtype).newbyteorder("<")
        if self.dtype not in CODE:
            raise TypeError("Prototype supports float64 or complex128")
        if not (0 <= lo <= hi <= global_size and 0 <= rank < ranks):
            raise ValueError("Invalid PETSc ownership range or MPI rank")
        if max_queued < 1:
            raise ValueError("max_queued must be >= 1")
        self.local_size = hi - lo
        self._socket = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        deadline = time.monotonic() + connect_timeout
        while True:
            try:
                self._socket.connect(str(socket_path))
                break
            except (FileNotFoundError, ConnectionRefusedError):
                if time.monotonic() >= deadline:
                    raise TimeoutError("Binary writer did not open socket")
                time.sleep(0.02)
        self._socket.sendall(HEADER.pack(
            MAGIC, rank, ranks, global_size, lo, hi, CODE[self.dtype]
        ))
        self._pending = queue.Queue(maxsize=max_queued)
        self._error = None
        self._closed = False
        self._thread = threading.Thread(
            target=self._drain, name=f"scrimp-trajectory-transport-{rank}",
            daemon=True,
        )
        self._thread.start()

    def _drain(self):
        try:
            while True:
                record = self._pending.get()
                try:
                    if record is None:
                        self._socket.sendall(FRAME.pack(END, 0, 0.))
                        return
                    step, t, data = record
                    self._socket.sendall(FRAME.pack(DATA, int(step), float(t)))
                    self._socket.sendall(memoryview(data).cast("B"))
                finally:
                    self._pending.task_done()
        except BaseException as exc:
            self._error = exc
        finally:
            self._socket.close()

    def enqueue(self, step, t, owned_state):
        """Copy only this rank's OWNED PETSc vector segment, without I/O."""
        if self._closed:
            raise RuntimeError("Stream is closed")
        if self._error is not None:
            raise RuntimeError("Binary writer transport failed") from self._error
        array = np.asarray(owned_state)
        if array.shape != (self.local_size,) or array.dtype != self.dtype:
            raise ValueError(
                "Rank data must be exact local PETSc array shape/dtype: "
                f"expected ({self.local_size},) {self.dtype}; "
                f"got {array.shape} {array.dtype}"
            )
        if not np.all(np.isfinite(array)):
            raise FloatingPointError("Cannot save nonfinite PETSc state")
        # Never retain a view onto PETSc TS-owned storage after callback.
        local_copy = np.array(array, copy=True, order="C")
        try:
            self._pending.put_nowait((int(step), float(t), local_copy))
        except queue.Full as exc:
            raise AsyncOutputBackpressure(
                "Asynchronous binary queue full: disk/transport cannot "
                "keep up with the simulation. No state was silently dropped. "
                "Increase max_queued, slow save cadence or use faster storage."
            ) from exc

    def close(self):
        """Drain all pending saves; waiting here is allowed after TS.solve."""
        if self._closed:
            return
        self._closed = True
        if self._error is not None:
            raise RuntimeError("Writer transport failed") from self._error
        self._pending.put(None)
        self._thread.join()
        if self._error is not None:
            raise RuntimeError("Writer transport failed") from self._error


class _Handler(socketserver.BaseRequestHandler):
    def handle(self):
        server = self.server
        try:
            header = HEADER.unpack(_recv_exact(self.request, HEADER.size))
            magic, rank, ranks, n_global, lo, hi, dtype_code = header
            if (magic != MAGIC or ranks != server.ranks
                    or n_global < 0 or not 0 <= rank < ranks
                    or not 0 <= lo <= hi <= n_global
                    or dtype_code not in DTYPES):
                raise ValueError(f"Invalid rank registration {header}")
            with server.lock:
                if rank in server.seen:
                    raise ValueError(f"Duplicate writer rank {rank}")
                server.seen.add(rank)
            dtype = DTYPES[dtype_code]
            nbytes = (hi - lo) * dtype.itemsize
            stem = f"trajectory_rank{rank:04d}.bin"
            partial = server.output_dir / (stem + ".partial")
            final = server.output_dir / stem
            samples = []
            last_step = None
            buf = bytearray(nbytes)
            with partial.open("wb") as output:
                output.write(HEADER.pack(*header))
                while True:
                    kind, step, t = FRAME.unpack(
                        _recv_exact(self.request, FRAME.size)
                    )
                    if kind == END:
                        break
                    if kind != DATA or not np.isfinite(t):
                        raise ValueError("Invalid trajectory frame")
                    if last_step is not None and step == last_step:
                        raise ValueError("Duplicate simulation step")
                    last_step = step
                    # Consume a complete frame before acknowledging its
                    # validity. The producer is decoupled through a bounded
                    # queue/thread and never blocks on these disk writes.
                    if nbytes:
                        buf[:] = _recv_exact(self.request, nbytes)
                    output.write(FRAME.pack(kind, step, t))
                    output.write(buf)
                    samples.append((step, t))
                output.flush()
                os.fsync(output.fileno())
            os.replace(partial, final)
            with server.lock:
                server.done[rank] = {
                    "rank": rank, "global_size": n_global, "lo": lo,
                    "hi": hi, "dtype_code": dtype_code,
                    "filename": stem, "samples": samples,
                }
        except BaseException as exc:
            with server.lock:
                server.errors.append(repr(exc))


class _Server(socketserver.ThreadingUnixStreamServer):
    daemon_threads = True
    allow_reuse_address = False

    def __init__(self, socket_path, output_dir, ranks):
        super().__init__(str(socket_path), _Handler)
        self.output_dir = Path(output_dir)
        self.ranks = ranks
        self.lock = threading.Lock()
        self.seen = set()
        self.done = {}
        self.errors = []


def serve(socket_path, output_dir, ranks, timeout=60.):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    socket_path = Path(socket_path)
    if socket_path.exists():
        raise FileExistsError(f"Refusing to reuse socket: {socket_path}")
    with _Server(socket_path, output_dir, ranks) as server:
        server.timeout = 0.2
        deadline = time.monotonic() + timeout
        worker = threading.Thread(target=server.serve_forever, daemon=True)
        worker.start()
        try:
            while True:
                with server.lock:
                    if server.errors:
                        raise RuntimeError(f"Writer failure: {server.errors}")
                    if len(server.done) == ranks:
                        completed = dict(server.done)
                        break
                if time.monotonic() >= deadline:
                    raise TimeoutError("Not all ranks finished writing")
                time.sleep(0.025)
        finally:
            server.shutdown()
            worker.join()
    try:
        socket_path.unlink()
    except FileNotFoundError:
        pass

    parts = sorted(completed.values(), key=lambda record: record["lo"])
    expected_lo = 0
    first = parts[0]
    reference_samples = first["samples"]
    for part in parts:
        if part["lo"] != expected_lo:
            raise ValueError("Rank ownership has gap or overlap")
        expected_lo = part["hi"]
        if (part["dtype_code"] != first["dtype_code"]
                or part["global_size"] != first["global_size"]
                or part["samples"] != reference_samples):
            raise ValueError("Ranks disagree on scalar type, size or samples")
    if expected_lo != first["global_size"]:
        raise ValueError("Rank ownership does not cover full state")

    manifest = {
        "format": "scrimp-rank-binary-v1",
        "ranks": ranks,
        "global_size": first["global_size"],
        "scalar_type": DTYPES[first["dtype_code"]].name,
        "number_of_states": len(reference_samples),
        "steps": [step for step, t in reference_samples],
        "times": [t for step, t in reference_samples],
        "partitions": [
            {key: part[key] for key in ("rank", "lo", "hi", "filename")}
            for part in parts
        ],
    }
    partial = output_dir / "trajectory_manifest.json.partial"
    final = output_dir / "trajectory_manifest.json"
    partial.write_text(json.dumps(manifest, indent=2) + "\n")
    os.replace(partial, final)
    return final


class LazyTrajectory:
    """Read one global sample at a time from per-rank memory-mapped files."""

    def __init__(self, manifest_path):
        manifest_path = Path(manifest_path)
        data = json.loads(manifest_path.read_text(encoding="utf-8"))
        if data["format"] != "scrimp-rank-binary-v1":
            raise ValueError("Unsupported trajectory format")
        self.global_size = data["global_size"]
        self.times = np.asarray(data["times"], dtype=np.float64)
        self.steps = np.asarray(data["steps"], dtype=np.int64)
        self.dtype = np.dtype(data["scalar_type"])
        if self.dtype not in CODE:
            raise ValueError("Unsupported stored dtype")
        self._parts = []
        for part in data["partitions"]:
            nlocal = part["hi"] - part["lo"]
            record_dtype = np.dtype([
                ("frame", "u1"), ("step", "<i8"),
                ("time", "<f8"), ("values", self.dtype, (nlocal,)),
            ], align=False)
            file = manifest_path.parent / part["filename"]
            # Framing on disk is FRAME size: one uint8 + int64 + float64.
            if record_dtype.itemsize != FRAME.size + nlocal * self.dtype.itemsize:
                raise AssertionError("Invalid binary record layout")
            mm = np.memmap(
                file, mode="r", dtype=record_dtype, offset=HEADER.size,
                shape=(len(self.times),),
            )
            self._parts.append((part["lo"], part["hi"], mm))

    def __len__(self):
        return len(self.times)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return np.stack([self[i] for i in range(*index.indices(len(self)))])
        index = int(index)
        if index < 0:
            index += len(self)
        if not (0 <= index < len(self)):
            raise IndexError(index)
        out = np.empty(self.global_size, dtype=self.dtype)
        for lo, hi, array in self._parts:
            out[lo:hi] = array["values"][index]
        return out


def start_writer(output_dir, ranks, *, timeout=60.):
    """Start a FRESH interpreter, never an MPI fork containing PETSc."""
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    socket_path = output_dir / "trajectory_writer.sock"
    logfile = (output_dir / "trajectory_writer.log").open("w")
    try:
        proc = subprocess.Popen(
            [sys.executable, str(Path(__file__).resolve()), "--writer",
             "--socket", str(socket_path), "--output-dir",
             str(output_dir), "--ranks", str(ranks),
             "--timeout", str(timeout)],
            stdin=subprocess.DEVNULL, stdout=logfile,
            stderr=subprocess.STDOUT, close_fds=True,
        )
    finally:
        logfile.close()
    return proc, socket_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--writer", action="store_true")
    parser.add_argument("--socket", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--ranks", type=int)
    parser.add_argument("--timeout", type=float, default=60.)
    args = parser.parse_args()
    if not args.writer or not args.socket or not args.output_dir or not args.ranks:
        parser.error("Writer mode needs --writer, --socket, --output-dir, --ranks")
    serve(args.socket, args.output_dir, args.ranks, args.timeout)


if __name__ == "__main__":
    main()
