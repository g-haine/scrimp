# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Compare saved *native* original shallow-water production physical TS MPI1/2.

No Docker, GetFEM, PETSc or MPI needed. Checks actual accepted BDF2
time steps, published initial nonlinear_affine algebraic residual,
full-state and Hamiltonian rank parity, and explicit gyro accepted-step
freeze contract. Does NOT assert exact energy/power preservation for
BDF2 with lagged explicit gyro.
"""

import argparse
import json
import math
from pathlib import Path


def require(value, message):
    if not value:
        raise ValueError(message)


def square_norm(sequence):
    return math.sqrt(math.fsum(float(x)*float(x) for x in sequence))


def check(folder):
    folder = Path(folder)
    reports = {}
    for ranks in (1, 2):
        f = folder/f"result_{ranks}"/f"production_affine_physical_{ranks}r.json"
        require(f.is_file(), f"Missing native JSON: {f}")
        data = json.loads(f.read_text(encoding="utf-8"))
        require(data.get("schema") == "scrimp-v130-issue28-production-affine-short-physical-ts-v1",
                "Wrong source evidence schema")
        require(data.get("status") == "PASS_NATIVE_PHYSICAL_TS_SMOKE",
                "Native physical-time result was not a PASS")
        require(data.get("mpi_ranks") == ranks, "Incorrect MPI rank count")
        require(data.get("source") == "original published examples/shallow_water.py",
                "Source physical model was modified")
        require(data.get("original_mesh") == .1, "Original mesh was changed")
        require(data.get("physical_initial_h0_p0_unchanged") is True,
                "Original differential initial data changed")
        require(data.get("original_homogeneous_boundary_controls_unchanged") is True,
                "Physical boundary control changed")
        require(data.get("initialization_strategy") == "nonlinear_affine",
                "Wrong production DAE initializer selected")
        init = data.get("production_initialization_report", {})
        require(init.get("strategy") == "nonlinear_affine"
                and init.get("status") == "PASS",
                "Real production opt-in initial projection not accepted")
        require(init.get("differential_initial_exact") is True,
                "Initial h0,p0 were modified")
        require(init.get("accuracy_target_met") is True,
                "Initial algebraic acceptance was not reached")
        require(init.get("algebraic_null_columns") == 0,
                "Structural singularity was overlooked")
        require(init.get("linear_solver_converged_reason", 0) > 0,
                "Initial direct solver did not converge")
        require(init.get("dofs") == 4270
                and init.get("differential_dofs") == 2035
                and init.get("algebraic_dofs") == 2235,
                "Original finite-element system structure changed")
        fa = float(init.get("native_final_algebraic_residual_l2", math.inf))
        require(math.isfinite(fa) and fa <= 1e-11,
                "Published algebraic residual exceeds strict target")
        physical = data["physical"]
        require(physical.get("manufactured_nonzero_G_initial_data") is False
                and physical.get("p_initial_expression") == "[0.,0.]",
                "Original physical initial p0 was altered")
        steps = physical.get("time_steps")
        require(steps == 4 and physical.get("TS_accepted_steps") == steps,
                "No four accepted true physical PETSc TS steps")
        require(physical.get("TS_converged_reason", 0) > 0
                and physical.get("TS_rejected_steps") == 0
                and physical.get("TS_SNES_failures") == 0,
                "Physical BDF2 solve failed or rejected steps")
        require(physical.get("time_scheme") == "bdf", "Physical BDF2 scheme changed")
        require(physical.get("max_frozen_explicit_RHS_drift") == 0,
                "Frozen accepted-step gyro changed during Newton")
        require(physical.get("main_IFunction_checks", 0) >= steps,
                "Physical IFunction not invoked")
        times = physical.get("times", [])
        Z = physical.get("full_state_trajectory", [])
        H = physical.get("hamiltonian", [])
        require(len(times) == len(Z) == len(H) >= 2,
                "Missing full PETSc physical-time trajectory")
        require(all(math.isfinite(x) for x in times + H),
                "Non-finite times or Hamiltonian")
        require(all(0 < h for h in H), "Hamiltonian nonpositive")
        require(all(x < y for x, y in zip(times, times[1:])),
                "Physical trajectory times are not strictly increasing")
        require(all(len(row) == 4270 and all(math.isfinite(v) for v in row)
                    for row in Z), "Nonfinite or truncated state samples")
        require(math.isclose(times[-1], steps*physical["dt"], rel_tol=0, abs_tol=1e-10),
                "Final physical time does not match accepted steps")
        events = physical.get("main_explicit_prestep", [])
        require(len(events) == steps, "Wrong accepted-step explicit cache events")
        require(all(math.isclose(e["time"], e["frozen_rhs_time"],
                     rel_tol=0, abs_tol=1e-12) for e in events),
                "Explicit cache not frozen at correct accepted time")
        require(data.get("exact_Hamiltonian_balance_qualified") is False,
                "Exact Hamiltonian power balance wrongly claimed")
        require(data.get("rotating_high_index_initialization_qualified") is False,
                "Rotating high-index Y0 wrongly marked qualified")
        reports[ranks] = data

    one = reports[1]["physical"]
    two = reports[2]["physical"]
    t1,t2 = one["times"],two["times"]
    require(len(t1) == len(t2), "Different MPI trajectory lengths")
    require(all(math.isclose(a,b,rel_tol=0,abs_tol=1e-12)
                for a,b in zip(t1,t2)), "MPI time samples differ")
    max_z = 0.
    max_h = 0.
    for a,b in zip(one["full_state_trajectory"],two["full_state_trajectory"]):
        err = square_norm([x-y for x,y in zip(a,b)])
        normalized = err/max(1.,square_norm(a),square_norm(b))
        max_z = max(max_z, normalized)
    for a,b in zip(one["hamiltonian"],two["hamiltonian"]):
        max_h=max(max_h, abs(a-b)/max(1.,abs(a),abs(b)))
    require(max_z <= 5e-8 and max_h <= 5e-8,
            f"MPI trajectory parity failed: state={max_z:.9g}, H={max_h:.9g}")
    result={
        "schema":"scrimp-v130-issue28-production-affine-physical-comparison-v1",
        "status":"PASS_PHYSICAL_TS_SMOKE_MPI_PARITY",
        "actual_petsc_physical_accepted_steps":4,
        "initialization_strategy":"nonlinear_affine",
        "original_h0_p0_and_boundary_controls_preserved":True,
        "mpi_max_scaled_full_state_difference":max_z,
        "mpi_max_relative_hamiltonian_difference":max_h,
        "signed_hamiltonian_change_mpi1":one["signed_H_final_minus_H_initial"],
        "signed_hamiltonian_change_mpi2":two["signed_H_final_minus_H_initial"],
        "exact_energy_power_balance_claimed":False,
        "rotating_high_index_DAE_initialization_claimed":False,
    }
    (folder/"production_affine_short_physical_comparison.json").write_text(
        json.dumps(result,indent=2)+"\n", encoding="utf-8"
    )
    print(
        "PASS: real production nonlinear_affine initial projection + 4 true "
        "BDF2 physical TS steps, native MPI 1/2 parity",
        flush=True,
    )
    print(
        f"full-state MPI relative difference={max_z:.9g}; "
        f"Hamiltonian MPI relative difference={max_h:.9g}; "
        f"measured signed H drift rank1={one['signed_H_final_minus_H_initial']:.9g}",
        flush=True,
    )
    print(
        "Hamiltonian power budget and rotating high-index traction Y0 "
        "are NOT yet qualified.",
        flush=True,
    )
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--evidence",type=Path,required=True)
    args=ap.parse_args()
    check(args.evidence)


if __name__=="__main__":
    main()
