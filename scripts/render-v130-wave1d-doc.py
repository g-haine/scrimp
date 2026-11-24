#!/usr/bin/env python3
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Render authentic 1D wave GetFEM VTK + CN evidence for Sphinx and users.

Input: the genuine five-solver native acceptance archive (ZIP or directory)
produced by scripts/run-v130-first-wave-1d-acceptance.sh.
Output: reproducible SVG/PNG figures, an animated GIF, useful CSV tables,
a compact NPZ and provenance JSON.

The animation's displacement w is RECONSTRUCTED: q=dw/dx is integrated
exactly on each P2 element; the left displacement is obtained by
trapezoidal time integration of the observed left velocity with gauge
w(0,0)=0. This is an illustration, NOT an independent conservation law.
Do not apply a scalar field offset or claim that q is displacement.
"""

import argparse
import base64
import csv
import hashlib
import io
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import zipfile

import numpy as np


VALID_CASES = (
    "mumps_mpi1", "mumps_mpi2", "superlu_mpi1",
    "superlu_dist_mpi1", "superlu_dist_mpi2",
)
NATIVE_PASS = "PASS_FIRST_WAVE_1D_ALL_AVAILABLE_LU_MPI_PHYSICS_VTK"


class Evidence:
    """Read a real on-disk acceptance archive, without extracting its files."""

    def __init__(self, path):
        self.path = Path(path)
        if self.path.is_file() and self.path.suffix.lower() == ".zip":
            self.archive = zipfile.ZipFile(self.path)
        elif self.path.is_dir():
            self.archive = None
        else:
            raise FileNotFoundError(path)

    def read(self, relative):
        relative = Path(relative)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("Unsafe evidence filename")
        if self.archive is not None:
            return self.archive.read(relative.as_posix())
        return (self.path / relative).read_bytes()


def vtk_float32(array):
    """Read the uncompressed GetFEM binary Float32 DataArray payload."""
    if (array is None or array.attrib.get("format") != "binary"
            or array.attrib.get("type") != "Float32"):
        raise ValueError("Expected binary GetFEM Float32 VTK DataArray")
    data = base64.b64decode("".join((array.text or "").split()), validate=True)
    if len(data) < 8 or (len(data)-4) % 4:
        raise ValueError("Bad GetFEM Float32 binary payload length")
    count = (len(data)-4) // 4
    length_bytes = int.from_bytes(data[:4], "little")
    float_count = float(np.frombuffer(data[:4], dtype="<f4")[0])
    if length_bytes != 4*count and float_count != count:
        raise ValueError("Invalid GetFEM binary VTK length header")
    return np.frombuffer(data[4:], dtype="<f4").astype(float)


def native_vtu(payload, variable):
    piece = ET.fromstring(payload).find(".//Piece")
    if piece is None:
        raise ValueError("Missing VTK unstructured grid piece")
    n = int(piece.attrib["NumberOfPoints"])
    ncells = int(piece.attrib["NumberOfCells"])
    if n < 2 or ncells < 1:
        raise ValueError("Incomplete VTK grid")
    coordinates = vtk_float32(piece.find("./Points/DataArray"))
    if len(coordinates) != 3*n:
        raise ValueError("Incorrect VTK coordinate array")
    points = coordinates.reshape(-1, 3)
    if np.max(np.abs(points[:, 1:])) > 1e-7:
        raise ValueError("Expected a 1D spatial mesh")
    values = vtk_float32(
        piece.find(f"./PointData/DataArray[@Name='{variable}']")
    )
    if len(values) != n or not np.all(np.isfinite(values)):
        raise ValueError("Missing or nonfinite original VTK field")
    order = np.argsort(points[:, 0])
    x, y = points[order, 0], values[order]
    if (not np.all(np.diff(x) > 0) or abs(x[0]) > 1e-6
            or abs(x[-1]-1.) > 1e-6):
        raise ValueError("Original 1D VTK mesh is incomplete/unsorted")
    return x, y, ncells


def read_run(data, case):
    """Validate five-solver acceptance before preparing visual figures."""
    if case not in VALID_CASES:
        raise ValueError("Not a qualified direct-solver case")
    summary = json.loads(data.read("summary.json"))
    if summary["status"] != NATIVE_PASS:
        raise ValueError("No valid complete native wave 1D acceptance")
    if summary["cases"][case]["status"] != "PASS_NATIVE_REAL_DAE_CN_ENERGY_AND_VTK":
        raise ValueError("Chosen native run was not accepted")
    with np.load(io.BytesIO(data.read(f"{case}/rank0.npz")),
                 allow_pickle=False) as bundle:
        trace = {name: bundle[name].copy() for name in bundle.files}
    times = trace["times"]
    if (times.ndim != 1 or len(times) < 2
            or not np.all(np.diff(times) > 0)):
        raise ValueError("Invalid accepted FE simulation times")
    original = {}
    for field in ("q", "p"):
        directory = f"{case}/pv/{field}"
        entries = ET.fromstring(
            data.read(f"{directory}/{field}.pvd")
        ).findall(".//DataSet")
        if len(entries) != len(times):
            raise ValueError(f"Original {field} PVD time series incomplete")
        grids, arrays = [], []
        for k, entry in enumerate(entries):
            if abs(float(entry.attrib["timestep"]) - times[k]) > 1e-10:
                raise ValueError("PVD time differs from accepted CN time")
            name = entry.attrib["file"]
            if name != Path(name).name:
                raise ValueError("Unsafe PVD VTU reference")
            x, values, _ = native_vtu(data.read(f"{directory}/{name}"), field)
            grids.append(x)
            arrays.append(values)
        if any(not np.array_equal(grids[0], v) for v in grids[1:]):
            raise ValueError("Changing GetFEM mesh during fixed simulation")
        original[field] = (grids[0], np.stack(arrays))
    xq, q = original["q"]
    xp, p = original["p"]
    if (len(xq) != 2*len(xp)-1
            or not np.allclose(xq[::2], xp, atol=1e-7, rtol=0)):
        raise ValueError("Expected compatible P2 strain / P1 momentum grids")

    # Physical velocity v=p/rho at the P2 output points. Interpolate
    # P1 momentum to xq BEFORE dividing by the analytical density.
    p_on_q = np.stack([np.interp(xq, xp, v) for v in p])
    velocity = p_on_q/(1.+xq*(1.-xq))[None, :]

    # Exact spatial integrals of piecewise quadratic Lagrange strain:
    # half interval: h*(5*q0+8*qm-q1)/24
    # full interval: h*(q0+4*qm+q1)/6.
    relative_w = np.zeros_like(q)
    for j in range(len(xp)-1):
        i = 2*j
        h = xq[i+2]-xq[i]
        relative_w[:, i+1] = (
            relative_w[:, i] + h*(5*q[:, i]+8*q[:, i+1]-q[:, i+2])/24.
        )
        relative_w[:, i+2] = (
            relative_w[:, i] + h*(q[:, i]+4*q[:, i+1]+q[:, i+2])/6.
        )
    left_velocity = p[:, 0]  # rho(0)=1.
    left_w = np.zeros(len(times))
    left_w[1:] = np.cumsum(
        np.diff(times)*(left_velocity[:-1]+left_velocity[1:])/2.
    )
    displacement = relative_w + left_w[:, None]
    fields = dict(
        t=times, xq=xq, xp=xp, q=q, p=p,
        velocity=velocity, w=displacement, w_left=left_w,
    )
    if any(not np.all(np.isfinite(v)) for v in fields.values()):
        raise ValueError("Nonfinite derived physical field")
    return summary, trace, fields


def save_png_svg(fig, output):
    fig.savefig(str(output)+".svg", bbox_inches="tight")
    fig.savefig(str(output)+".png", dpi=160, bbox_inches="tight")
    import matplotlib.pyplot as plt
    plt.close(fig)


def make_plots(out, trace, fields, summary, case):
    import matplotlib.pyplot as plt

    t, x = fields["t"], fields["xq"]
    q, v, w = fields["q"], fields["velocity"], fields["w"]
    chosen = np.unique(np.linspace(0, len(t)-1, min(5, len(t))).astype(int))

    fig, axs = plt.subplots(3, 1, figsize=(9, 9), sharex=True,
                            layout="constrained")
    for k in chosen:
        for ax, var in zip(axs, (q, v, w)):
            ax.plot(x, var[k], label=f"t={t[k]:.3f} s", lw=1.8)
    for ax, title in zip(
        axs, ("Strain q = dx w (P2)",
              "Velocity v = p/rho (P1 momentum)",
              "Reconstructed displacement w (gauge w(0,0)=0)"),
    ):
        ax.set_title(title, loc="left")
        ax.grid(True, alpha=.3)
    axs[0].legend(ncol=min(5, len(chosen)), fontsize=8)
    axs[-1].set_xlabel("Position x in [0,1]")
    fig.suptitle(f"String: actual GetFEM field snapshots ({case})")
    save_png_svg(fig, out/"wave_1d_field_profiles")

    fig, axs = plt.subplots(3, 1, figsize=(9, 8), sharex=True,
                            layout="constrained")
    for ax, var, name in zip(
        axs, (q, v, w), ("Strain q", "Velocity p/rho",
                       "Reconstructed displacement w"),
    ):
        plot = ax.pcolormesh(x, t, var, shading="auto")
        fig.colorbar(plot, ax=ax, shrink=.8)
        ax.set_ylabel("Time (s)")
        ax.set_title(name, loc="left")
    axs[-1].set_xlabel("Position x")
    fig.suptitle("Wave propagation: accepted FE states through time")
    save_png_svg(fig, out/"wave_1d_spacetime")

    H = trace["H"]
    W = trace["signed_work"]
    fig, axs = plt.subplots(3, 1, figsize=(9, 8), sharex=True,
                            layout="constrained")
    axs[0].plot(t, H-H[0], "o-", label="H(t)-H(0)")
    axs[0].plot(t, W, "--", label="CN-integrated boundary work")
    axs[0].legend()
    axs[0].set_ylabel("Energy and work")
    axs[1].plot(t, trace["power_left"], label="Left boundary power")
    axs[1].plot(t, trace["power_right"], label="Right boundary power")
    axs[1].legend()
    axs[1].set_ylabel("Signed power")
    axs[2].plot(t, H-H[0]-W, "o-")
    axs[2].set_ylabel("H-H0-work")
    axs[2].ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    axs[2].set_xlabel("Time (s)")
    for ax in axs:
        ax.grid(True, alpha=.3)
    fig.suptitle("Mechanical power balance verified with real CN powers")
    save_png_svg(fig, out/"wave_1d_energy_balance")

    parity = summary["comparisons"]
    names = [x for x in parity if x != "mumps_mpi1"]
    fig, ax = plt.subplots(figsize=(9, 4), layout="constrained")
    xx = np.arange(len(names))
    for j, (key, caption) in enumerate((
        ("states", "Full state"), ("H", "Hamiltonian"),
        ("power_left", "Left port power"),
    )):
        y = [max(float(parity[name][key]), 1e-18) for name in names]
        ax.bar(xx+(j-1)*.24, y, width=.24, label=caption)
    ax.set_xticks(xx, [x.replace("_", " ") for x in names])
    ax.set_yscale("log")
    ax.set_ylabel("Max abs difference vs MUMPS, MPI1")
    ax.legend()
    ax.grid(True, axis="y", alpha=.25)
    ax.set_title("Five solvers/MPI configurations on one identical physical model")
    save_png_svg(fig, out/"wave_1d_solver_parity")


def animate(out, fields, case, milliseconds=250):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    t, x = fields["t"], fields["xq"]
    values = (fields["q"], fields["velocity"], fields["w"])
    fig, axs = plt.subplots(3, 1, figsize=(8, 7), sharex=True,
                            layout="constrained")
    lines = []
    for ax, arr, caption in zip(
        axs, values, ("Strain q", "Velocity p/rho",
                      "Reconstructed displacement w"),
    ):
        line, = ax.plot(x, arr[0], lw=2)
        lines.append(line)
        lo, hi = float(arr.min()), float(arr.max())
        margin = max(.025, .08*(hi-lo))
        ax.set_ylim(lo-margin, hi+margin)
        ax.set_ylabel(caption)
        ax.grid(True, alpha=.25)
    axs[-1].set_xlabel("Position x")
    title = fig.suptitle("")

    def draw(k):
        for line, arr in zip(lines, values):
            line.set_ydata(arr[k])
        title.set_text(f"SCRIMP 1D string • {case} • t={t[k]:.3f} s")
        return (*lines, title)

    movie = FuncAnimation(fig, draw, frames=len(t),
                          interval=milliseconds, blit=False)
    movie.save(out/"wave_1d_string.gif",
               writer=PillowWriter(fps=1000./milliseconds), dpi=110)
    plt.close(fig)


def export_data(out, trace, fields, manifest, case, source_sha256):
    np.savez_compressed(
        out/"wave_1d_postprocessed.npz", **fields,
        H=trace["H"], power_left=trace["power_left"],
        power_right=trace["power_right"],
        signed_work=trace["signed_work"],
        energy_defect=trace["energy_defect"],
    )
    with (out/"wave_1d_energy_and_ports.csv").open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow((
            "time_s", "H", "H_minus_H0", "signed_boundary_work",
            "power_left", "power_right", "H_minus_H0_minus_work",
        ))
        for k, t in enumerate(fields["t"]):
            writer.writerow((
                float(t), float(trace["H"][k]),
                float(trace["H"][k]-trace["H"][0]),
                float(trace["signed_work"][k]),
                float(trace["power_left"][k]),
                float(trace["power_right"][k]),
                float(trace["H"][k]-trace["H"][0]-trace["signed_work"][k]),
            ))
    x, xp = fields["xq"], fields["xp"]
    with (out/"wave_1d_fields.csv").open("w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow((
            "time_s", "x", "strain_q", "momentum_p_interpolated",
            "velocity_p_over_rho", "displacement_w_reconstructed",
        ))
        for k, t in enumerate(fields["t"]):
            pi = np.interp(x, xp, fields["p"][k])
            for j, coordinate in enumerate(x):
                writer.writerow((
                    float(t), float(coordinate), float(fields["q"][k, j]),
                    float(pi[j]), float(fields["velocity"][k, j]),
                    float(fields["w"][k, j]),
                ))
    metadata = {
        "status": "RENDERED_FROM_GENUINE_ACCEPTED_NATIVE_FE_AND_VTK",
        "source_commit": manifest["source"],
        "qualified_image_revision": manifest["image_revision"],
        "input_SHA256": source_sha256,
        "native_case": case,
        "frames": len(fields["t"]),
        "P2_nodes": len(fields["xq"]),
        "P1_nodes": len(fields["xp"]),
        "max_original_CN_step_energy_defect":
            float(np.max(np.abs(trace["energy_defect"]))),
        "displacement_is_reconstructed": True,
        "initial_displacement_gauge": "w(0,0)=0",
        "spatial_quadrature": "exact P2 polynomial element integration",
        "left_boundary_motion": "trapezoidal integration of p(0,t)/rho(0)",
        "precision_warning":
            "VTU PointData is Float32; visualisation only. "
            "Original Hamiltonian and port powers come from real "
            "GetFEM FE measurements, not the approximate visual grid.",
        "not_an_independent_physical_qualification": True,
    }
    (out/"provenance.json").write_text(
        json.dumps(metadata, indent=2)+"\n", encoding="utf-8"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True,
                        help="Native accepted five-backend ZIP or its directory")
    parser.add_argument("--output", type=Path, required=True,
                        help="Host-owned docs / user image output directory")
    parser.add_argument("--case", choices=VALID_CASES, default="mumps_mpi1")
    parser.add_argument("--frame-ms", type=int, default=250)
    parser.add_argument("--no-gif", action="store_true")
    args = parser.parse_args()
    if not 70 <= args.frame_ms <= 3000:
        parser.error("Animation frame duration must be 70..3000 ms")
    data = Evidence(args.input)
    summary, trace, fields = read_run(data, args.case)
    manifest = json.loads(data.read("manifest.json"))
    args.output.mkdir(parents=True, exist_ok=True)
    sha = (hashlib.sha256(args.input.read_bytes()).hexdigest()
           if args.input.is_file() else None)
    export_data(args.output, trace, fields, manifest, args.case, sha)
    make_plots(args.output, trace, fields, summary, args.case)
    if not args.no_gif:
        animate(args.output, fields, args.case, args.frame_ms)
    print(json.dumps({
        "status": "RENDERED_REAL_NATIVE_WAVE_1D",
        "source_commit": manifest["source"],
        "case": args.case,
        "frames": len(fields["t"]),
        "output": str(args.output),
        "assets": sorted(p.name for p in args.output.iterdir()),
    }, indent=2))


if __name__ == "__main__":
    main()
