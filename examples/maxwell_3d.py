# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
"""Small genuine 3D Maxwell cavity using patched intrinsic GetFEM elements.

Normalized lossless Maxwell: eps=mu=1, no currents or port drive.
E (electric field) in patched Nedelec-II k=1 H(curl), B (magnetic
induction) in patched tetrahedral BDM1 H(div), in a real 3D ball.
The original native Nedelec0/RT0 baseline at commit daac18cb
is preserved as its separately qualified archived simulation.

The weak equations are
    (dE/dt, v) = (B, curl v),
    (dB/dt, w) = -(curl E, w).
The first equation has natural zero tangential magnetic boundary work;
this is a weak homogeneous closure, not a separately enforced PEC
tangential E condition. The discrete coupling is skew-adjoint, hence
the *exact semidiscrete and CN identity* H(t)=H(0), with
H=0.5 integral(|E|^2+|B|^2).

Run only inside the pinned SCRIMP GetFEM 5.5 MPI image:
    python -m examples.maxwell_3d --check --export --output /artifacts
This example is not an electromagnetic convergence or divergence
benchmark. It deliberately tests the Nedelec and RT GetFEM patches.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import scrimp as S

def assert_patched_intrinsic_fem(fem, expected_type, expected_order):
    """Reject the exact legacy S.FEM fallback without MeshFem.fem() introspection.

    SCRIMP's FEM.set_fem() explicitly changes its *own* type and order
    to CG/1 if native GetFEM rejects the requested family. We retain
    the same FEM object that was passed to DPHS.add_FEM(), so its
    post-installation metadata is authoritative for that fallback.
    Its real GetFEM MeshFem must also have qdim=3 and nonzero DOFs.
    """
    actual_type = fem.get_type()
    actual_order = fem.get_order()
    if actual_type != expected_type or actual_order != expected_order:
        raise RuntimeError(
            f"GetFEM intrinsic {expected_type}/{expected_order} not active: "
            f"SCRIMP changed to {actual_type}/{actual_order}"
        )
    mesh_fem = fem.get_fem()
    if mesh_fem is None:
        raise RuntimeError(f"Missing GetFEM MeshFem for {expected_type}")
    try:
        qdim = int(mesh_fem.qdim())
        ndof = int(mesh_fem.nbdof())
    except (TypeError, AttributeError, ValueError) as exc:
        raise RuntimeError(
            f"Cannot verify GetFEM MeshFem dimension for {expected_type}"
        ) from exc
    if qdim != 3 or ndof <= 0:
        raise RuntimeError(
            f"Invalid intrinsic {expected_type} 3D field: qdim={qdim}, DOFs={ndof}"
        )
    return ndof


# TUTORIAL M3D-01: explicit three-dimensional curl (GWFL has no native Curl)
CURL3 = (
    "[Grad(v)(3,2)-Grad(v)(2,3);"
    "Grad(v)(1,3)-Grad(v)(3,1);"
    "Grad(v)(2,1)-Grad(v)(1,2)]"
)




def maxwell_3d(*, mesh_size=0.7, final_time=0.2, time_step=0.02,
               factor_solver="mumps"):
    """Build and solve the two-field Maxwell system on a tetrahedral ball."""
    values = (mesh_size, final_time, time_step)
    if any(not np.isfinite(float(v)) or float(v) <= 0 for v in values):
        raise ValueError("Mesh and time parameters must be finite and positive")
    if abs(final_time/time_step-round(final_time/time_step)) > 1e-10:
        raise ValueError("Choose an integer number of fixed CN steps")
    if factor_solver not in ("mumps", "superlu", "superlu_dist"):
        raise ValueError("Unknown PETSc LU backend")

    model = S.DPHS("real")
    model.gf_model.add_macro("Curl3(v)", CURL3)

    # TUTORIAL M3D-02: real 3D tetrahedral mesh, provided by SCRIMP
    model.set_domain(S.Domain("Ball", {"R": 1., "h": float(mesh_size)}))
    if model.domain.get_dim()[0] != 3:
        raise RuntimeError("Maxwell requires a genuine 3D volume mesh")

    # TUTORIAL M3D-03: coenergy = state (eps=mu=1)
    E = S.State("E", "Electric field", "vector-field")
    B = S.State("B", "Magnetic induction", "vector-field")
    model.add_state(E)
    model.add_state(B)
    model.add_costate(S.CoState("e_E", "Electric field", E, substituted=True))
    model.add_costate(S.CoState("e_B", "Magnetic induction", B, substituted=True))

    # TUTORIAL M3D-04: exercise BOTH GetFEM source patches, not built-in
    # Nedelec0/RT0. NED2 k=1 -> FEM_NEDELEC2(3,1) (patch 0004),
    # BDM k=1 -> FEM_BDMK(3,1) (tetrahedral DoF patch 0001).
    # curl([P1]^3) is constant and is contained in BDM1;
    # exact div(curl(E)) = 0 is preserved by the mass equation.
    fem_E = S.FEM("E", 1, "NED2")
    fem_B = S.FEM("B", 1, "BDM")

    # S.FEM has a legacy silent CG fallback on native AssertionError.
    # Retain these exact objects to detect it after registration.
    for fem, expected_type, expected_order in (
        (fem_E, "NED2", 1),
        (fem_B, "BDM", 1),
    ):
        model.add_FEM(fem)
        assert_patched_intrinsic_fem(fem, expected_type, expected_order)

    # TUTORIAL M3D-05: two positive mass blocks, two skew curl blocks
    for brick in (
        S.Brick("M_E", "E.Test_E", [1], dt=True, position="flow"),
        S.Brick("M_B", "B.Test_B", [1], dt=True, position="flow"),
        S.Brick("C^T", "B.Curl3(Test_E)", [1], position="effort"),
        S.Brick("-C", "-Curl3(E).Test_B", [1], position="effort"),
    ):
        model.add_brick(brick)

    # TUTORIAL M3D-06: smooth solenoidal rotation (exact Nedelec-II P1 mode),
    # zero magnetic field; curl E=[0,0,2] triggers nontrivial dynamics.
    model.set_initial_value("E", "[-y, x, 0.]")
    model.set_initial_value("B", "[0., 0., 0.]")

    model.set_time_scheme(
        # Both states carry a positive mass: this is an ODE, NOT a DAE.
        # The algebraic initializer requires algebraic unknowns;
        # a pseudo-TS or artificial DAE spin-up would be wrong here.
        init_step=False,
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        t_0=0., t_f=float(final_time),
        dt=float(time_step), dt_save=float(time_step),
        ksp_type="gmres", pc_type="lu",
        pc_factor_mat_solver_type=factor_solver,
    )
    model.solve()

    # TUTORIAL M3D-07: true GetFEM 3D electromagnetic Hamiltonian
    model.hamiltonian.set_name("3D electromagnetic energy")
    model.hamiltonian.add_term(S.Term("Electric", "0.5*E.E", [1]))
    model.hamiltonian.add_term(S.Term("Magnetic", "0.5*B.B", [1]))
    return model


def measure_maxwell(system):
    """Audit real FE energy, nontrivial B, and discrete Gauss diagnostic."""
    import getfem as gf

    t = np.asarray(system.solution["t"], dtype=float)
    z = np.asarray(system.solution["z"], dtype=float)
    H = np.asarray(system.get_Hamiltonian(), dtype=float)
    if t.ndim != 1 or len(t) < 2 or z.ndim != 2 or len(z) != len(t):
        raise RuntimeError("Invalid original FE Maxwell time/state arrays")
    if H.shape != t.shape or not np.all(np.diff(t) > 0):
        raise RuntimeError("Invalid Maxwell energy sampling")
    if not all(np.all(np.isfinite(x)) for x in (t, z, H)):
        raise RuntimeError("Nonfinite Maxwell FE data")
    mim = system.domain._int_method[0]
    electric, magnetic, gauss_B = [], [], []
    for tk, zk in zip(t, z):
        system.gf_model.set_time(float(tk))
        system.gf_model.to_variables(zk)
        electric.append(float(gf.asm(
            "generic", mim, 0, "0.5*E.E", 1, system.gf_model)))
        magnetic.append(float(gf.asm(
            "generic", mim, 0, "0.5*B.B", 1, system.gf_model)))
        gauss_B.append(float(gf.asm(
            "generic", mim, 0,
            "Trace(Grad(B))*Trace(Grad(B))", 1, system.gf_model)))
    He = np.asarray(electric)
    Hb = np.asarray(magnetic)
    div_B_squared = np.asarray(gauss_B)
    steps = np.diff(H)
    defect = H-H[0]
    scale = max(1.0, float(H[0]))
    gates = {
        "real_3d_mesh": system.domain.get_dim()[0] == 3,
        "authentic_fixed_CN_steps": bool(
            system.last_ts_converged_reason > 0
            and system.last_ts_steps == len(t)-1
            and system.last_ts_rejections == 0
            and system.last_ts_snes_failures == 0),
        "pure_differential_Maxwell_no_artificial_DAE": bool(
            len(system.states) == 2 and
            all(not p.get_algebraic() for p in system.ports.values())),
        "positive_initial_em_energy": bool(
            H[0] > 1e-7 and min(np.min(He), np.min(Hb)) >= -1e-10),
        "electric_and_magnetic_exchange": bool(
            Hb[-1] > 1e-6*H[0] and He[-1] < He[0]),
        "independent_fe_energy_matches": bool(
            np.max(np.abs(He+Hb-H)) < 1e-9*scale),
        "every_CN_step_conserves_energy": bool(
            np.max(np.abs(steps)) < 1e-9*scale),
        "cumulative_em_energy_conserved": bool(
            np.max(np.abs(defect)) < 1e-9*scale),
        # Gauss(B): curl E belongs to RT0 on the same tetrahedral mesh;
        # B(0)=0, hence the compatible discrete div B must stay zero.
        # This guards the *second* Maxwell constraint independently
        # of the skew/energy check, not only as a plotted diagnostic.
        "magnetic_gauss_law": bool(
            np.min(div_B_squared) >= -1e-12*scale and
            np.max(div_B_squared) <= 1e-8*scale),
    }
    return dict(t=t, z=z, H=H, electric_energy=He, magnetic_energy=Hb,
                div_B_squared=div_B_squared, step_energy_defect=steps,
                cumulative_energy_defect=defect, gates=gates)



# Visualization is deliberately kept IN the canonical executable.
# Only rank 0 calls these routines, AFTER collective GetFEM VTK export.
# The electric/magnetic energy curves use original GetFEM FE integrals;
# all 3-D arrows and the wireframe use native GetFEM VTK snapshots.
def _vtk_native_array(element, expected_count=None):
    """Read GetFEM's inline-base64 VTK, without requiring the vtk package."""
    import base64
    import struct

    if element is None or element.attrib.get("format") != "binary":
        raise ValueError("Missing binary GetFEM VTK DataArray")
    dtype = {"Float32": "<f4", "Float64": "<f8",
             "Int32": "<i4", "Int64": "<i8"}.get(element.get("type"))
    if dtype is None:
        raise ValueError("Unsupported GetFEM VTK DataArray type")
    raw = base64.b64decode((element.text or "").strip(), validate=True)
    if len(raw) < 4:
        raise ValueError("Truncated GetFEM VTK binary array")
    if expected_count is None:
        if struct.unpack("<I", raw[:4])[0] != len(raw) - 4:
            raise ValueError("Invalid GetFEM VTK binary byte-size header")
    elif len(raw) != 4 + expected_count * np.dtype(dtype).itemsize:
        # PointData's 4-byte header can encode a component count rather
        # than a byte length. Trust the Piece's explicit field dimensions.
        raise ValueError("GetFEM VTK field length disagrees with mesh")
    a = np.frombuffer(raw[4:], dtype=np.dtype(dtype)).copy()
    if expected_count is not None and len(a) != expected_count:
        raise ValueError("VTK numerical field has incorrect size")
    return a


def _native_vector_history(out, field, times):
    """Read ONLY the real exported tetrahedral vector snapshots."""
    from xml.etree import ElementTree as ET

    folder = Path(out) / "pv" / field
    catalog = ET.parse(folder / f"{field}.pvd").getroot()
    frames = catalog.findall(".//DataSet")
    if len(frames) != len(times):
        raise ValueError(f"{field}: VTK frame count differs from FE trajectory")
    points0, tetra0, snapshots = None, None, []
    for i, frame in enumerate(frames):
        if abs(float(frame.attrib["timestep"]) - float(times[i])) > 1e-10:
            raise ValueError(f"{field}: VTK/FE time mismatch at frame {i}")
        filename = Path(frame.attrib["file"])
        if filename.name != str(filename) or filename.suffix != ".vtu":
            raise ValueError(f"{field}: unsafe or invalid VTK frame filename")
        tree = ET.parse(folder / filename).getroot()
        piece = tree.find(".//Piece")
        if piece is None:
            raise ValueError(f"{field}: missing VTK unstructured grid")
        npnt = int(piece.attrib["NumberOfPoints"])
        nc = int(piece.attrib["NumberOfCells"])
        if npnt <= 0 or nc <= 0:
            raise ValueError(f"{field}: empty native 3D VTK mesh")
        points = _vtk_native_array(
            piece.find("./Points/DataArray"), expected_count=npnt*3
        ).reshape(npnt, 3)
        conn = _vtk_native_array(
            piece.find("./Cells/DataArray[@Name='connectivity']")
        )
        offsets = _vtk_native_array(
            piece.find("./Cells/DataArray[@Name='offsets']")
        )
        types = _vtk_native_array(
            piece.find("./Cells/DataArray[@Name='types']")
        )
        if (len(conn) != 4*nc or not np.array_equal(
            offsets, 4*np.arange(1, nc+1)
        ) or len(types) != nc or not np.all(types == 10)):
            raise ValueError(f"{field}: non-tetrahedral or corrupt VTK cells")
        tetra = conn.reshape(nc, 4).astype(int)
        if np.any((tetra < 0) | (tetra >= npnt)):
            raise ValueError(f"{field}: invalid tetrahedron vertex index")
        node = piece.find(f"./PointData/DataArray[@Name='{field}']")
        if node is None or int(node.get("NumberOfComponents", "0")) != 3:
            raise ValueError(f"{field}: expected native 3-vector field")
        values = _vtk_native_array(node, expected_count=npnt*3).reshape(
            npnt, 3
        )
        if not np.all(np.isfinite(points)) or not np.all(np.isfinite(values)):
            raise ValueError(f"{field}: nonfinite GetFEM VTK data")
        if points0 is None:
            points0, tetra0 = points, tetra
        elif not (np.array_equal(points, points0)
                  and np.array_equal(tetra, tetra0)):
            raise ValueError(f"{field}: VTK mesh geometry changed with time")
        snapshots.append(values)
    return points0, tetra0, np.asarray(snapshots)


def render_maxwell_visuals(output, observed, state_sha256):
    """Create PNG/SVG/GIF INSIDE this runnable example, from real FE/VTK.

    Called only by MPI rank 0 after both native fields are exported.
    It does not change GetFEM state, timestepping or the physical gates.
    """
    import collections
    import io
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from PIL import Image

    output = Path(output)
    folder = output / "visuals"
    folder.mkdir(parents=True, exist_ok=True)
    time = np.asarray(observed["t"], dtype=float)
    E_pts, E_cells, E_vals = _native_vector_history(output, "E", time)
    B_pts, B_cells, B_vals = _native_vector_history(output, "B", time)
    if not (np.array_equal(E_pts, B_pts)
            and np.array_equal(E_cells, B_cells)):
        raise ValueError("Electric/magnetic VTK exports do not share a mesh")

    # Only edges from the *outer* tetrahedral surface are drawn.
    faces = collections.Counter(
        tuple(sorted(face)) for a, b, c, d in E_cells
        for face in ((a,b,c), (a,b,d), (a,c,d), (b,c,d))
    )
    edges = sorted({
        tuple(sorted((u, v)))
        for face, count in faces.items() if count == 1
        for u, v in ((face[0],face[1]), (face[0],face[2]),
                     (face[1],face[2]))
    })
    if not edges:
        raise ValueError("Cannot identify exterior tetrahedral surface")
    radius = max(1., float(np.max(np.abs(E_pts))) * 1.1)

    fig, ax = plt.subplots(figsize=(8.8, 5.0))
    ax.plot(time, observed["electric_energy"], marker="o",
            label="Electric energy")
    ax.plot(time, observed["magnetic_energy"], marker="s",
            label="Magnetic energy")
    ax.plot(time, observed["H"], marker=".", label="Hamiltonian")
    ax.set(xlabel="Time (dimensionless)", ylabel="Energy (dimensionless)",
           title="Maxwell 3D: electric–magnetic energy exchange\n"
                 "Nédélec II (k=1) / BDM1 · native GetFEM quadratures")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(folder / "maxwell3d_echange_energie.png", dpi=170)
    fig.savefig(folder / "maxwell3d_echange_energie.svg")
    plt.close(fig)

    def vector_figure(field, vectors, frame):
        fig = plt.figure(figsize=(7.0, 6.6))
        ax = fig.add_subplot(111, projection="3d", proj_type="ortho")
        for i, j in edges:
            ax.plot(E_pts[[i,j],0], E_pts[[i,j],1], E_pts[[i,j],2],
                    linewidth=0.5, alpha=0.12)
        # At most 160 true VTK nodes for legibility; no resampling.
        selected = np.linspace(
            0, len(E_pts)-1, min(len(E_pts), 160), dtype=int
        )
        v = vectors[frame, selected]
        p = E_pts[selected]
        ax.quiver(p[:,0], p[:,1], p[:,2],
                  v[:,0], v[:,1], v[:,2],
                  length=0.45, normalize=False, arrow_length_ratio=0.22)
        ax.set(xlim=(-radius,radius), ylim=(-radius,radius),
               zlim=(-radius,radius), xlabel="x", ylabel="y", zlabel="z")
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=22, azim=-58)
        kind = "electric field E" if field == "E" else "magnetic induction B"
        ax.set_title(f"Maxwell 3D: {kind}\n"
                     f"t = {time[frame]:.2f} · native GetFEM 3D vectors",
                     pad=12)
        fig.subplots_adjust(left=0.03, right=0.97,
                            bottom=0.02, top=0.88)
        return fig

    for field, values in (("E",E_vals), ("B",B_vals)):
        fig = vector_figure(field, values, -1)
        fig.savefig(folder / f"maxwell3d_vecteurs_{field}_tfinal.png",
                    dpi=155)
        plt.close(fig)

    # Keep the same camera, domain bounds and vector scale in *all*
    # frames: an initially zero B field must really appear from zero.
    gif_frames = []
    for i in range(len(time)):
        fig = vector_figure("B", B_vals, i)
        buffer = io.BytesIO()
        fig.savefig(buffer, format="png", dpi=105)
        plt.close(fig)
        buffer.seek(0)
        with Image.open(buffer) as im:
            gif_frames.append(im.convert("RGB").copy())
    gif_frames[0].save(
        folder / "maxwell3d_B_3D_animation.gif",
        save_all=True, append_images=gif_frames[1:],
        duration=300, loop=0, optimize=False,
    )
    assets = (
        "maxwell3d_echange_energie.png",
        "maxwell3d_echange_energie.svg",
        "maxwell3d_vecteurs_E_tfinal.png",
        "maxwell3d_vecteurs_B_tfinal.png",
        "maxwell3d_B_3D_animation.gif",
    )
    provenance = {
        "schema": "scrimp-maxwell3d-native-visuals-v1",
        "origin": "examples/maxwell_3d.py (same run as FE solution)",
        "source": "native GetFEM 3D VTU/PVD and FE Hamiltonian arrays",
        "MPI_writer_rank": 0,
        "original_FE_state_SHA256": state_sha256,
        "frames": int(len(time)),
        "points": int(len(E_pts)),
        "tetrahedra": int(len(E_cells)),
        "fields": ["E", "B"],
        "no_synthetic_field_or_time_interpolation": True,
        "visual_note": "VTK nodal samples represent underlying intrinsic FE; "
                       "arrows are not edge/face DOFs.",
        "assets": list(assets),
    }
    (folder / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    (folder / "README.txt").write_text(
        "Generated by examples/maxwell_3d.py from this very run.\n"
        "Energy: original GetFEM FE quadratures; arrows: original VTK nodal "
        "3D vector samples on the real tetrahedral mesh.\n"
        "The GIF uses fixed camera, vector scale and coordinates.\n"
        "Arrows are nodal visualizations, not Nedelec edge or BDM face DOFs.\n",
        encoding="utf-8",
    )
    return provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh-size", type=float, default=0.7)
    parser.add_argument("--time-step", type=float, default=0.02)
    parser.add_argument("--final-time", type=float, default=0.2)
    parser.add_argument("--solver", choices=("mumps", "superlu", "superlu_dist"),
                        default="mumps")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--export", action="store_true")
    args = parser.parse_args()
    if (args.check or args.export) and args.output is None:
        parser.error("--check/--export require --output")
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    if comm.size > 1 and args.solver == "superlu":
        parser.error("Sequential SuperLU cannot run in MPI2")

    system = maxwell_3d(mesh_size=args.mesh_size,
                        final_time=args.final_time,
                        time_step=args.time_step, factor_solver=args.solver)
    if args.output is None:
        return system
    args.output.mkdir(parents=True, exist_ok=True)
    obs = measure_maxwell(system)
    np.savez_compressed(
        args.output/f"rank{comm.rank}.npz",
        **{k: v for k, v in obs.items() if isinstance(v, np.ndarray)}
    )
    metadata = {
        "status": ("PASS_NATIVE_MAXWELL3D" if all(obs["gates"].values())
                   else "FAIL_NATIVE_MAXWELL3D"),
        "rank": comm.rank, "MPI": comm.size,
        "solver": args.solver,
        "FE": {
            "E": {"family": "NED2", "k": 1,
                  "getfem": "FEM_NEDELEC2(3,1)",
                  "DOFs": int(system.ports["E"].get_fem().nbdof())},
            "B": {"family": "BDM", "k": 1,
                  "getfem": "FEM_BDMK(3,1)",
                  "DOFs": int(system.ports["B"].get_fem().nbdof())},
        },
        "DOFs": int(obs["z"].shape[1]),
        "steps": int(system.last_ts_steps),
        "TS_reason": int(system.last_ts_converged_reason),
        "DAE_initialization": system.initialization_report,
        "physical_gates": obs["gates"],
        "initial_energy": float(obs["H"][0]),
        "final_energy": float(obs["H"][-1]),
        "final_magnetic_energy": float(obs["magnetic_energy"][-1]),
        "max_step_energy_defect": float(np.max(np.abs(
            obs["step_energy_defect"]))),
        "max_cumulative_energy_defect": float(np.max(np.abs(
            obs["cumulative_energy_defect"]))),
        "max_divB_squared": float(np.max(obs["div_B_squared"])),
        "state_SHA256": hashlib.sha256(
            np.ascontiguousarray(obs["z"]).tobytes()).hexdigest(),
    }
    (args.output/f"rank{comm.rank}.json").write_text(
        json.dumps(metadata, indent=2)+"\n", encoding="utf-8"
    )
    if args.check and not all(obs["gates"].values()):
        raise RuntimeError("Native Maxwell 3D physical acceptance FAILED: "
                           + json.dumps(metadata))
    if args.export:
        # Native VTK writer is collective; it serializes only on rank 0.
        for name in ("E", "B"):
            system.export_to_pv(name, path=str(args.output), t="All")
        # Figures/animation are local file operations, rank 0 ONLY.
        # --export therefore produces the visual package automatically.
        if comm.rank == 0:
            metadata["visuals"] = render_maxwell_visuals(
                args.output, obs, metadata["state_SHA256"]
            )
            (args.output/"rank0.json").write_text(
                json.dumps(metadata, indent=2)+"\n", encoding="utf-8"
            )
    if comm.rank == 0:
        print(json.dumps(metadata, indent=2), flush=True)
    return system


if __name__ == "__main__":
    main()
