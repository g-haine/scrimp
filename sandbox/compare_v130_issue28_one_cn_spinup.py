# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Host comparator for native #28/#6 one PETSc TSCN spinup step.

PETSc TSCN uses theta=0.5 ENDPOINT trapezoidal, NOT implicit midpoint.
Validate its actual equation F(t1,z1,2*delta/dt)+F(t0,z0,0)=0,
nonlinear and algebraic blocks, real TS acceptance and MPI full vectors.
The independent midpoint norm is OUT-OF-SCHEME information only.
PETSc warns that TSCN endpoint assumes constant linear mass whereas
the SWE mass is nonlinear M(h); time order/power remain UNQUALIFIED.
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np

OCI_REV = "1d6d0c90115213298047b28237a7f38087563a35"
MESHES = (.2, .1)


def require(predicate, message):
    if not predicate:
        raise ValueError(message)


def finite(value, name):
    require(type(value) in (int, float) and math.isfinite(value),
            f"Invalid nonfinite {name}: {value!r}")
    return float(value)


def validate(report, mesh, ranks):
    require(report.get("schema") == "scrimp-v130-issue28-one-cn-spinup-v1"
            and report.get("status") ==
            "NATIVE_SINGLE_PHYSICAL_CN_STEP_DIAGNOSTIC_ONLY",
            "Wrong native single-step diagnostic schema")
    require(report.get("source_mesh") == mesh
            and report.get("mpi_ranks") == ranks
            and report.get("new_physical_rest_spinup_experiment") is True
            and report.get("original_incompatible_rotation_initial_state_unchanged") is True
            and report.get("geometry_order") == 2
            and report.get("height_order") == 2
            and report.get("physical_facet_tangent_wall") is True
            and report.get("rotation_after_spinup") == .2
            and report.get("spinup_duration") == .02,
            "Original physical P2 or rest-spinup source changed")
    initial = report.get("initial", {})
    require(initial.get("strategy") ==
            "sandbox_exact_stationary_spinup_only"
            and initial.get("status") ==
            "SANDBOX_ONLY_PUBLISHED_NOT_GENERAL_INITIALIZER"
            and initial.get("differential_q0_preserved_bitwise") is True
            and initial.get("getfem_globalvec_published_bitwise") is True
            and initial.get("no_pseudo_or_BDF_initialization") is True
            and initial.get("analytic_gprime_t0_exact_zero") is True,
            "Unsafe/non-analytic initial-state publication")
    init_f = initial.get("initial_native_IFunction", {})
    require(finite(init_f.get("algebraic_l2"), "initial Fa") < 1e-10
            and finite(init_f.get("dynamic_l2"), "initial Fd") < 1e-10,
            "Physical t0 was not truly algebraically/dynamically consistent")
    require(report.get("physical_scheme") == "cn"
            and report.get("physical_time_step") == .002
            and report.get("ts_accepted_steps") == 1
            and report.get("ts_rejected_steps") == 0
            and report.get("ts_snes_failures") == 0
            and type(report.get("ts_converged_reason")) is int
            and report["ts_converged_reason"] > 0
            and abs(finite(report.get("ts_actual_final_time"), "TS time")-.002)
            < 1e-12,
            "Real PETSc CN did not accept exactly one physical step")
    require(report.get("two_endpoint_states_archived") is True
            and report.get("native_original_TSCN_endpoint_equation_verified") is True
            and report.get("implicit_midpoint_scheme_used") is False
            and report.get("petsc_TSCN_form") ==
            "theta_0p5_endpoint_trapezoidal"
            and report.get("petsc_TSCN_shift") == 1000.
            and report.get("petsc_TSCN_affine_old_F_coefficient") == -1.
            and report.get("nonlinear_mass_PETSc_endpoint_assumption_violated") is True
            and report.get("nonlinear_mass_time_accuracy_qualified") is False
            and report.get("production_high_index_initializer_qualified") is False
            and report.get("energy_balance_convergence_qualified") is False
            and report.get("physical_multistep_trajectory_qualified") is False,
            "Wrong PETSc TSCN method or unjustified nonlinear mass/energy qualification")
    cn = report.get("original_IFunction_TSCN_endpoint_equation_l2", {})
    mid = report.get("original_IFunction_midpoint_OUT_OF_SCHEME_l2", {})
    end = report.get("original_IFunction_endpoint_zero_dot_l2", {})
    original = report.get("original_IFunction_t0_zero_dot_l2", {})
    final_scaled = report.get("original_IFunction_t1_twice_difference_dot_l2", {})
    for key in ("full_l2", "dynamic_l2", "algebraic_l2"):
        require(finite(cn.get(key), "TSCN endpoint equation "+key) < 1e-7,
                "Accepted real TSCN does not satisfy exact PETSc endpoint equation")
    require(finite(end.get("algebraic_l2"), "endpoint algebraic") < 1e-7,
            "Accepted TSCN step violates endpoint algebraic constraint")
    for name, stats in (
        ("old zero-dot IFunction", original),
        ("new twice-difference-dot IFunction", final_scaled),
        ("out-of-scheme implicit midpoint", mid),
        ("new endpoint zero-dot IFunction", end),
    ):
        for key in ("algebraic_l2", "dynamic_l2", "full_l2"):
            require(finite(stats.get(key), name+" "+key) >= 0,
                    "Missing independent original IFunction evaluation")
    for key in ("U_l2",):
        require(finite(end.get(key), "endpoint "+key) >= 0,
                "Nonfinite original endpoint boundary constraint")

    for field in (
        "hamiltonian_initial", "hamiltonian_final",
        "hamiltonian_change", "physical_boundary_normal_work",
        "physical_boundary_viscous_work", "integrated_boundary_work",
        "integrated_viscous_dissipation", "integrated_dissipation_grad",
        "integrated_dissipation_div", "physical_power_balance_defect",
        "physical_power_balance_defect_over_actual_exchange",
        "physical_volume_balance_defect", "dynamic_step_change_l2",
        "first_step_minimum_height_coefficient",
    ):
        finite(report.get(field), field)
    require(report["first_step_minimum_height_coefficient"] > 0,
            "Physical h coefficients nonpositive")
    require(report.get("first_physical_state_changed") is True
            and report["dynamic_step_change_l2"] > 1e-12,
            "Accepted physical TS did not actually change dynamic h/p")

    require(report["integrated_viscous_dissipation"] >= -1e-12,
            "Physical viscous dissipation negative")
    require(math.isclose(
        report["integrated_boundary_work"],
        report["physical_boundary_normal_work"]
        + report["physical_boundary_viscous_work"],
        rel_tol=1e-7, abs_tol=1e-13,
    ), "Physical SWE wall power does not match weighted port work")
    require(math.isclose(
        report["integrated_viscous_dissipation"],
        report["integrated_dissipation_grad"]
        + report["integrated_dissipation_div"],
        rel_tol=1e-7, abs_tol=1e-13,
    ), "SWE viscous dissipation components do not sum")
    dh = report["hamiltonian_change"]
    dw = report["integrated_boundary_work"]
    diss = report["integrated_viscous_dissipation"]
    err = dh + diss - dw
    require(math.isclose(
        err, report["physical_power_balance_defect"],
        rel_tol=1e-7, abs_tol=1e-12,
    ), "Physical balance defect cannot be reproduced from saved terms")
    exchanged = max(abs(dh), abs(dw), abs(diss), 1e-14)
    require(math.isclose(
        err/exchanged,
        report["physical_power_balance_defect_over_actual_exchange"],
        rel_tol=1e-7, abs_tol=1e-8,
    ), "Actual exchanged-energy normalisation was misreported")
    return {
        "initial_Fa_l2": init_f["algebraic_l2"],
        "TSCN_endpoint_equation_l2": cn["full_l2"],
        "midpoint_OUT_OF_SCHEME_l2": mid["full_l2"],
        "endpoint_Fa_l2": end["algebraic_l2"],
        "endpoint_U_l2": end["U_l2"],
        "physical_power_defect": err,
        "power_defect_over_actual_exchange": err/exchanged,
        "physical_volume_defect": report["physical_volume_balance_defect"],
        "accepted_steps": 1,
    }


def compare(directory):
    root = Path(directory)
    manifest = json.loads((root/"manifest.json").read_text(encoding="utf-8"))
    require(manifest.get("schema") == "v130-issue28-one-cn-spinup-v1"
            and manifest.get("image_revision") == OCI_REV
            and manifest.get("physical_spinup_only") is True
            and manifest.get("production_source_unchanged") is True
            and manifest.get("one_physical_step_expected") is True,
            "Invalid single-step native gate provenance")
    summary = {}
    for mesh in MESHES:
        by_rank = []
        for ranks in (1, 2):
            f = root/f"result_h{str(mesh).replace('.', 'p')}_mpi{ranks}"/"one_cn_spinup.json"
            require(f.is_file(), f"Missing native one-CN evidence {f}")
            by_rank.append(validate(
                json.loads(f.read_text(encoding="utf-8")), mesh, ranks,
            ))
        a, b = by_rank
        for key in (
            "initial_Fa_l2", "TSCN_endpoint_equation_l2", "endpoint_Fa_l2",
            "endpoint_U_l2", "physical_power_defect",
            "physical_volume_defect",
        ):
            va, vb = finite(a[key], key), finite(b[key], key)
            # Norms can be ~1e-16; a blanket absolute tolerance 1e-9
            # would never detect a corrupt rank on the audited scale.
            floor = 5e-13 if "power" in key or "volume" in key else 5e-15
            require(abs(va-vb) <= floor + 1e-5*max(abs(va), abs(vb)),
                    f"Native MPI parity mismatch in {key}")
        endpoints = []
        for ranks in (1, 2):
            filename = root/(
                f"result_h{str(mesh).replace('.', 'p')}_mpi{ranks}"
            )/"one_cn_spinup_endpoints.npz"
            require(filename.is_file(),
                    f"Missing original native CN endpoints {filename}")
            with np.load(filename, allow_pickle=False) as z:
                require(set(z.files) >= {
                    "times", "initial", "final",
                    "H", "p_normal", "p_viscous", "d_grad", "d_div"
                }, "Missing archived original physical one-step states/powers")
                require(z["times"].shape == (2,)
                        and abs(float(z["times"][0])) < 1e-12
                        and abs(float(z["times"][1])-.002) < 1e-12,
                        "Invalid original physical two endpoint times")
                for key in ("initial", "final"):
                    require(z[key].ndim == 1
                            and z[key].size > 0
                            and np.all(np.isfinite(z[key])),
                            f"Invalid original physical {key} coefficients")
                endpoints.append({
                    key: np.asarray(z[key], dtype=float).copy()
                    for key in ("initial", "final")
                })
        for key in ("initial", "final"):
            x1, x2 = endpoints[0][key], endpoints[1][key]
            require(x1.shape == x2.shape,
                    "Native CN MPI endpoint dimensions changed")
            distance = float(np.linalg.norm(x1-x2))
            step_scale = float(np.linalg.norm(
                endpoints[0]["final"]-endpoints[0]["initial"]
            ))
            require(distance <= 1e-10+1e-5*step_scale,
                    f"Native CN MPI full {key} coefficient vector parity failed")
        require(np.array_equal(
            endpoints[0]["initial"], endpoints[1]["initial"]
        ), "MPI original published initial physical coefficients differ")
        key = f"h{mesh:g}"
        summary[key] = {
            "native_mpi1": a, "native_mpi2": b,
            "native_mpi_parity": True,
            "energy_balance_verified_with_dt_refinement": False,
        }
        print(
            f"h={mesh:g}: 1 accepted physical CN step MPI1/2 PASS "
            f"TSCN_eq={a['TSCN_endpoint_equation_l2']:.9g} "
            f"mid_OUT_OF_SCHEME={a['midpoint_OUT_OF_SCHEME_l2']:.9g} "
            f"end_Fa={a['endpoint_Fa_l2']:.9g} "
            f"energy_exchange_defect={a['power_defect_over_actual_exchange']:.9g}",
            flush=True
        )
    result = {
        "schema": "scrimp-v130-issue28-one-cn-comparison-v1",
        "status": "ONE_NATIVE_TSCN_STEP_ACCEPTED_NONLINEAR_MASS_AND_POWER_UNQUALIFIED",
        "cases": summary,
        "physical_energy_balance_refinement_qualified": False,
        "nontrivial_rotating_swe_trajectory_qualified": False,
        "general_high_index_initializer_qualified": False,
        "interpretation": (
            "One PETSc TSCN endpoint-trapezoidal step uses a sandbox "
            "analytic initial state. The actual original IFunction "
            "endpoint equation F(t1,z1,2*delta/dt)+F(t0,z0,0) is checked. "
            "The implicit midpoint is OUT-OF-SCHEME only. PETSc documents "
            "that this endpoint formula assumes linear constant mass, "
            "which is NOT valid for nonlinear M(h) of SWE. No physical "
            "time-order, power balance or general initializer is qualified."
        ),
    }
    (root/"one_cn_spinup_comparison.json").write_text(
        json.dumps(result, indent=2)+"\n", encoding="utf-8",
    )
    print("COMPLETE: PETSc TSCN equation checked; nonlinear-mass physical order and energy NOT qualified.")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--evidence", type=Path, required=True)
    args = p.parse_args()
    compare(args.evidence)


if __name__ == "__main__":
    main()
