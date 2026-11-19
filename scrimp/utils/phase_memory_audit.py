# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Opt-in, local-only resident-memory sampling at solver lifecycle phases.

No PETSc, GetFEM, MPI collectives, garbage collections, matrix copies or
forced malloc_trim are performed. These samples show *where* live RSS
changes, not the exact ownership of native allocations. Linux-only.
"""
import json
import os
from pathlib import Path
import resource
import time


_FIELDS = (
    "VmRSS", "VmHWM", "VmSize", "VmPeak", "RssAnon",
    "RssFile", "RssShmem", "VmSwap",
)


class PhaseMemoryAudit:
    def __init__(self, directory, *, rank, ranks, source_commit, mode):
        self.directory = Path(directory)
        self.rank = int(rank)
        self.ranks = int(ranks)
        self.source_commit = str(source_commit)
        self.mode = str(mode)
        self.samples = []
        self.start = time.perf_counter()

    def sample(self, phase, **metadata):
        stats = {}
        for line in Path("/proc/self/status").read_text(
            encoding="ascii"
        ).splitlines():
            label, _, raw = line.partition(":")
            if label not in _FIELDS:
                continue
            tokens = raw.split()
            if len(tokens) == 2 and tokens[1] == "kB":
                stats[label + "_kib"] = int(tokens[0])
        if "VmRSS_kib" not in stats or "VmHWM_kib" not in stats:
            raise RuntimeError("Missing Linux RSS/HWM process counters")
        self.samples.append({
            "phase": str(phase),
            "wall_since_start_s": time.perf_counter() - self.start,
            "ru_maxrss_kib": int(resource.getrusage(
                resource.RUSAGE_SELF
            ).ru_maxrss),
            **{key: float(val) if isinstance(val, float)
               else val for key, val in metadata.items()},
            **stats,
        })

    def write(self):
        self.directory.mkdir(parents=True, exist_ok=True)
        target = self.directory / f"memory_phase_rank{self.rank:04d}.json"
        scratch = target.with_suffix(".json.partial")
        obj = {
            "schema": "scrimp-phase-memory-v1",
            "rank": self.rank, "ranks": self.ranks,
            "source_commit": self.source_commit,
            "matrix_memory_mode": self.mode,
            "samples": self.samples,
        }
        with scratch.open("w", encoding="utf-8") as handle:
            json.dump(obj, handle, indent=2, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(scratch, target)
        return target
