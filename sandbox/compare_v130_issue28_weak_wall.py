# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Host comparison of native #28 boundary weak moments and mass Riesz norms.

Never treat a nonzero strong L2 wall mismatch, or an apparent small
moment, as proof of full nonlinear DAE initial compatibility. The
original rotating h0,p0, wall control, rest control, P2 and MPI1/2
must be preserved. Only read archived GetFEM native diagnostics.
"""

import argparse
import json
import math
from pathlib import Path

EXPECTED_IMAGE = "1d6d0c90115213298047b28237a7f38087563a35"
COMPONENTS = ("full", "normal", "tangential")


def require(ok, message):
    if not ok:
        raise ValueError(message)


def number(value, label):
    require(type(value) in (float, int) and math.isfinite(value),
            f"Invalid real {label}: {value!r}")
    return float(value)


def l2(vec):
    return math.sqrt(math.fsum(number(v, "vector")**2 for v in vec))


def difference(a, b):
    require(len(a) == len(b), "Native FE wall vector size mismatch")
    return l2([x-y for x, y in zip(a, b)])


def check_case(report, case, mesh, ranks):
    require(report.get("schema") == "scrimp-v130-issue28-weak-wall-moments-v1"
            and report.get("status") == "NATIVE_WEAK_MOMENT_DIAGNOSTIC_ONLY",
            "Wrong native wall moment report")
    require(report.get("mpi_ranks") == ranks
            and abs(number(report.get("nominal_mesh"), "h") - mesh) < 1e-13
            and report.get("rotation") == (0. if case == "rest" else .2),
            f"Wrong wall experiment {mesh}/{case}/{ranks}")
    require(report.get("geometry_order") == 2
            and report.get("original_h0_p0_Y_preserved_bitwise") is True
            and report.get("entire_state_rollback_bitwise") is True
            and report.get("no_initialized_traction_Y0") is True
            and report.get("no_initialization_or_physical_TS") is True
            and report.get("algebraic_feasibility_proven") is False,
            "Unqualified physical initialization or altered original problem")
    count = report.get("boundary_U_dofs")
    require(type(count) is int and count > 0, "Invalid U trace DOFs")

    raw_perimeters = report["raw_per_rank_perimeter"]
    require(len(raw_perimeters) == ranks, "Missing MPI perimeter contributions")
    raw_squared = report["raw_per_rank_strong_wall_squared_integral"]
    require(len(raw_squared) == ranks, "Missing strong mismatch contributions")
    p = number(report.get("boundary_perimeter"), "perimeter")
    strong = number(report.get("strong_wall_vector_mismatch_rms"), "strong RMS")
    require(p > 0 and strong >= 0, "Invalid physical L2 measures")
    require(report["perimeter_reduction_mode"] in (
        "serial", "replicated-global", "summed-local"
    ), "Invalid scalar MPI perimeter mode")
    require(report["strong_squared_integral_reduction_mode"] in (
        "serial", "replicated-global", "summed-local"
    ), "Invalid scalar MPI mismatch mode")
    for label, arr, saved, kind in (
        ("perimeter", raw_perimeters, p, report["perimeter_reduction_mode"]),
        ("squared mismatch", raw_squared, p*strong**2,
         report["strong_squared_integral_reduction_mode"]),
    ):
        values = [number(x, label) for x in arr]
        recomputed = values[0] if kind == "replicated-global" else sum(values)
        require(math.isclose(recomputed, saved, rel_tol=2e-9,
                             abs_tol=1e-12),
                f"Wrong GetFEM scalar MPI {label} reconstruction")
    # Verify full=normal+tangential moments in the *actual* U test space.
    moments = report["moment_vectors"]
    riesz = report["boundary_mass_riesz"]
    coefficient_vectors = {}
    for component in COMPONENTS:
        e = moments[component]
        v = e.get("assembled_vector")
        require(isinstance(v, list) and len(v) == count,
                f"Missing native {component} U moment coefficients")
        v = [number(x, f"{component} moment") for x in v]
        coefficient_vectors[component] = v
        raw_vectors = e.get("per_rank_native_vector")
        require(isinstance(raw_vectors, list) and len(raw_vectors) == ranks,
                f"Missing {component} GetFEM MPI vector contributions")
        each = []
        for raw in raw_vectors:
            require(len(raw) == count, "Wrong native U vector length")
            each.append([number(x, component) for x in raw])
        mode = e.get("mpi_vector_reduction_mode")
        require(mode in ("serial", "replicated-global", "summed-local"),
                "Wrong native weak vector assembly mode")
        if mode == "replicated-global":
            require(all(
                difference(v, each[0]) <=
                max(1e-12, 1e-7*l2(each[0]))
                for v in each[1:]
            ), f"Replicated GetFEM MPI weak vector mismatch: {component}")
        computed = each[0] if mode == "replicated-global" else [
            sum(x) for x in zip(*each)
        ]
        require(difference(v, computed) <= 1e-12 * max(1., l2(v)),
                f"GetFEM weak boundary MPI {component} reduction mismatch")
        require(math.isclose(
            l2(v), number(e.get("coefficient_l2"), "coefficient norm"),
            rel_tol=1e-8, abs_tol=1e-12,
        ), f"Wrong native {component} coefficient norm")
        wr = riesz[component]
        dual = number(wr.get("riesz_boundary_RMS"), "weak dual RMS")
        projected = number(
            wr.get("riesz_L2_projected_mismatch_l2"), "projected L2"
        )
        require(dual >= 0 and projected >= 0
                and math.isclose(projected, dual*math.sqrt(p),
                                 rel_tol=1e-8, abs_tol=1e-12),
                "Incorrect Riesz mass physical norm")
        require(number(
            wr.get("boundary_mass_solver_relative_residual"),
            "boundary mass residual",
        ) < 1e-8, "Boundary Riesz mass solve failed")
        require(wr.get("boundary_mass_is_negative_native_J_UU") is True,
                "Native boundary mass not used")

    mismatch = difference(
        coefficient_vectors["full"], [
            x+y for x, y in zip(coefficient_vectors["normal"],
                               coefficient_vectors["tangential"])
        ]
    )
    require(mismatch < 1e-10 * max(1., l2(coefficient_vectors["full"])),
            "Weak normal/tangential boundary moment decomposition fails")
    require(math.isclose(
        mismatch,
        number(report.get("normal_tangential_decomposition_l2"),
               "decomposition error"),
        rel_tol=1e-4, abs_tol=1e-12,
    ), "Archived weak decomposition error false")
    require(riesz["full"]["riesz_boundary_RMS"] <= strong + 1e-8,
            "Riesz L2 projection exceeds strong physical wall mismatch")

    if case == "rest":
        require(strong < 1e-12, "Exact rest strongly mismatched")
        require(l2(coefficient_vectors["full"]) < 1e-11
                and riesz["full"]["riesz_boundary_RMS"] < 1e-11,
                "Exact rest weak wall control failed")
    return {
        "rms_strong": strong,
        "rms_weak_riesz": riesz["full"]["riesz_boundary_RMS"],
        "rms_weak_normal": riesz["normal"]["riesz_boundary_RMS"],
        "rms_weak_tangent": riesz["tangential"]["riesz_boundary_RMS"],
        "coefficients": coefficient_vectors,
        "perimeter": p,
        "vector_reduction_mode":
            {component: moments[component]["mpi_vector_reduction_mode"]
             for component in COMPONENTS},
    }


def compare(evidence):
    folder = Path(evidence)
    manifest = json.loads((folder/"manifest.json").read_text(encoding="utf-8"))
    require(manifest.get("schema") == "v130-issue28-weak-wall-moments-v1"
            and manifest.get("image_revision") == EXPECTED_IMAGE,
            "Wrong wall moment campaign/image")
    require(manifest.get("original_fixed_h0_p0") is True
            and manifest.get("physical_initial_Y0_accepted") is False
            and manifest.get("no_physical_time_steps") is True,
            "Wall gate claims false initial state acceptance")

    summary = {}
    for mesh in (.2, .1):
        for case in ("rest", "rotating"):
            metrics = []
            for ranks in (1, 2):
                p = folder / (
                    f"result_h{str(mesh).replace('.', 'p')}_{case}_{ranks}"
                ) / "weak_wall.json"
                require(p.is_file(), f"Missing native report {p}")
                r = json.loads(p.read_text(encoding="utf-8"))
                metrics.append(check_case(r, case, mesh, ranks))
            a, b = metrics
            require(math.isclose(a["perimeter"], b["perimeter"],
                                 rel_tol=1e-9, abs_tol=1e-10),
                    "GetFEM MPI scalar perimeter parity failed")
            for label in ("rms_strong", "rms_weak_riesz",
                          "rms_weak_normal", "rms_weak_tangent"):
                require(math.isclose(a[label], b[label],
                                     rel_tol=1e-5, abs_tol=1e-10),
                        f"MPI physical weak/strong norm parity failed: {label}")
            for component in COMPONENTS:
                vec_a = a["coefficients"][component]
                vec_b = b["coefficients"][component]
                require(difference(vec_a, vec_b) <=
                        max(1e-10, 1e-5*l2(vec_a)),
                        f"MPI original weak vector parity failed: {component}")
            key = f"h{mesh:g}_{case}"
            summary[key] = {
                "original_strong_wall_mismatch_RMS": a["rms_strong"],
                "original_weak_wall_Riesz_RMS": a["rms_weak_riesz"],
                "normal_weak_Riesz_RMS": a["rms_weak_normal"],
                "tangential_weak_Riesz_RMS": a["rms_weak_tangent"],
                "getfem_parallel_vector_modes": b["vector_reduction_mode"],
                "mpi_native_weak_moment_parity": True,
                "weak_moments_numeric_nonzero": bool(
                    a["rms_weak_riesz"] > 1e-11
                ),
            }
            print(
                f"h={mesh:g} {case}: direct weak Riesz wall RMS="
                f"{a['rms_weak_riesz']:.9g} "
                f"(strong={a['rms_strong']:.9g}); "
                "native MPI1/2 parity PASS", flush=True
            )
    result = {
        "schema": "scrimp-v130-issue28-weak-wall-summary-v1",
        "status": "COMPLETE_WEAK_WALL_DIAGNOSTIC_NOT_INITIALIZER",
        "cases": summary,
        "rotating_exact_weak_DAE_compatibility_certified": False,
        "physical_initial_Y0_accepted": False,
        "physical_TS_executed": False,
        "interpretation": (
            "Native GetFEM first-order moments of the original physical "
            "wall mismatch, calibrated against the actual port-U Riesz "
            "mass. Weak mismatch needs an independent exact FE "
            "elimination argument before classifying the original "
            "high-index constrained DAE initial state."
        ),
    }
    (folder/"weak_wall_comparison.json").write_text(
        json.dumps(result, indent=2)+"\n", encoding="utf-8"
    )
    print("COMPLETE: original weak wall measured; rotating DAE/Y0 NOT initialized.")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, required=True)
    args = parser.parse_args()
    compare(args.evidence)


if __name__ == "__main__":
    main()
