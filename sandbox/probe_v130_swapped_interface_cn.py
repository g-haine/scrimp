# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Native strict-DAE/CN energy audit for Heat=annulus, Wave=disk.

Runs the CURRENT examples/heat_wave.py, including the corrected outward
interface laws, with the opt-in automatic interface compiler.
No physical sign patches, manual face regions, or boundary-data patches.
In contrast to the heat=1/wave=2 case, the exterior Dirichlet boundary
belongs to HEAT; do NOT clamp the disk wave velocity artificially.

An initial-data/DAE failure is preserved as evidence, not recast as a
pass. True scientific qualification requires an actual completed CN
trajectory, nonzero interface exchange, and small separate balance
defects.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from examples.heat_wave import heat_wave_eq
from measure_v130_hw_dual_facet_physics import (
    collect_heat_wave_dual_facet_physics,
)


def _numpy_vec(x):
    a = np.asarray(x, dtype=float)
    if not np.all(np.isfinite(a)):
        raise RuntimeError("Nonfinite physical trajectory")
    return a


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if (MPI.COMM_WORLD.size != 1 or PETSc.COMM_WORLD.getSize() != 1
            or os.environ.get("SCRIMP_V130_AUTO_INTERFACE_PILOT") != "1"):
        raise RuntimeError("Requires automatic interface compiler with MPI1")
    if args.out.exists() and any(args.out.iterdir()):
        raise RuntimeError("Refuse to overwrite nonempty evidence directory")
    args.out.mkdir(parents=True, exist_ok=True)

    # The real example builds every GetFEM FEM/brick, projects all three
    # raw initial fields, initializes the algebraic pH-DAE, executes
    # four genuine fixed CN steps, then plots and exports the fields.
    system = heat_wave_eq(
        heat_region=2, wave_region=1, factor_solver="mumps",
        final_time=0.004, time_step=0.001, save_step=0.001,
    )
    compiler = system._internal_interface_compiler
    if compiler is None or not compiler.completed:
        raise RuntimeError("Original swapped model did not use auto compiler")
    manifest = compiler.report
    if len(manifest["interfaces"]) != 1:
        raise RuntimeError("Expected one automatically oriented interface")
    side = manifest["interfaces"][0]
    if (side["logical_region"] != 10
            or side["subdomains"] != [1, 2]
            or side["facet_count"] != 78
            or side["port_owners"] != {
                "Interface Heat": 2, "Interface Wave": 1,
            }):
        raise RuntimeError("Wrong swapped geometric/interface ownership")

    heat = system.ports["Interface Heat"]
    wave = system.ports["Interface Wave"]
    heat_face = int(heat.get_assembly_region())
    wave_face = int(wave.get_assembly_region())
    if (heat.get_region() != 10 or wave.get_region() != 10
            or heat_face == wave_face
            or heat_face != side["side_regions"][2]
            or wave_face != side["side_regions"][1]):
        raise RuntimeError("One-sided orientation or logical label corrupted")

    # Original 0-input external boundary is now thermally owned.
    boundary = system.ports["Boundary"]
    if boundary.get_region() != 20 or boundary.get_position() != "effort":
        raise RuntimeError("Swapped exterior heat Dirichlet causality changed")

    dae = system.initialization_report
    if (not isinstance(dae, dict)
            or dae.get("status") != "PASS"
            or not dae.get("differential_initial_exact")
            or dae.get("pseudo_ts_executed")):
        raise RuntimeError("Strict algebraic initial projection did not pass")

    times = _numpy_vec(system.solution["t"]).reshape(-1)
    states = _numpy_vec(system.solution["z"])
    H = _numpy_vec(system.get_Hamiltonian()).reshape(-1)
    if (times.shape != (5,) or not np.allclose(
            times, np.arange(5) * 0.001, rtol=0, atol=1e-12)
            or states.shape != (5, 8056) or H.shape != (5,)
            or system.last_ts_converged_reason <= 0
            or system.last_ts_steps != 4
            or system.last_ts_rejections != 0
            or system.last_ts_snes_failures != 0):
        raise RuntimeError("Actual swapped fixed CN trajectory invalid")

    # Unlike the original 1/2 collector, all volumetric diagnostic
    # integrals must now be taken on heat=2 / wave=1.
    traces, physics = collect_heat_wave_dual_facet_physics(
        system, times, H,
        heat_faces=heat_face, wave_faces=wave_face,
        heat_bulk=2, wave_bulk=1,
    )
    power_checks = {}
    for port, expected_key in (
        (heat, "CN_P_interface_heat"),
        (wave, "CN_P_interface_wave"),
    ):
        if port.get_power_sampling() != "cn_step":
            raise RuntimeError("Public interface port not sampled at CN midpoints")
        actual = _numpy_vec(port.get_power()).reshape(-1)
        expected = _numpy_vec(traces[expected_key]).reshape(-1)
        if actual.shape != expected.shape:
            raise RuntimeError("Public interface port-power length mismatch")
        err = float(np.max(np.abs(actual - expected)))
        if not np.allclose(actual, expected, rtol=1e-10, atol=1e-10):
            raise RuntimeError("Public port power differs from oriented trace")
        power_checks[port.get_name()] = err

    boundary_power = float(physics["max_outer_port_power"])
    scale = max(1., abs(float(H[0])))
    power_scale = max(
        1., float(physics["max_CN_heat_dissipation"]),
        float(physics["max_individual_CN_interface_power"]),
    )
    gates = {
        "strict_initial_DAE": True,
        "four_steps_no_reject": True,
        "total_CN_balance":
            physics["max_CN_step_H_balance_defect"] <= 1e-7 * scale,
        "heat_CN_balance":
            physics["max_CN_step_heat_balance_defect"] <= 1e-7 * scale,
        "wave_CN_balance":
            physics["max_CN_step_wave_balance_defect"] <= 1e-7 * scale,
        "closed_interface_power":
            physics["max_CN_interface_power_defect"] <= 1e-7 * power_scale,
        "nonzero_heat_dissipation":
            physics["max_CN_heat_dissipation"] > 1e-12,
        "nontrivial_interface_excitation":
            physics["interface_exchange_nontrivial"]
            and physics["max_individual_CN_interface_power"] > 1e-8,
        "homogeneous_exterior_power":
            boundary_power <= 1e-8 * power_scale,
        "nonincreasing_H":
            physics["max_CN_positive_H_increment"] <= 1e-7 * scale,
    }
    passed = all(gates.values())
    result = {
        "schema": "scrimp-v130-swapped-auto-HW-real-CN-physics-v1",
        "status": (
            "PASS_SWAPPED_AUTO_INTERFACE_STRICT_DAE_CN_ENERGY_MPI1"
            if passed else
            "CN_COMPLETED_SWAPPED_PHYSICAL_BALANCE_NOT_QUALIFIED"
        ),
        "heat_bulk": 2,
        "wave_bulk": 1,
        "heat_faces": heat_face,
        "wave_faces": wave_face,
        "outer_thermal_boundary": 20,
        "initial_fields": "original L2-projected Gaussian and zero stress",
        "boundary_compatibility_patch_applied": False,
        "original_live_example_source_modified_by_test": False,
        "no_manual_interface_region_or_transmission_sign_override": True,
        "actual_4_step_CN": True,
        "full_state_shape": list(states.shape),
        "strict_DAE_report": dae,
        "geometric_compiler": manifest,
        "physical_gates": gates,
        "physics": physics,
        "public_Port_compute_CN_max_difference": power_checks,
        "MPI2_qualified": False,
    }
    (args.out / "swapped_auto_cn_audit.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8",
    )
    np.savez_compressed(
        args.out / "swapped_auto_cn_trajectory.npz",
        times=times, states=states, H=H,
        **{"phys_" + name: values for name, values in traces.items()},
    )
    print(json.dumps({
        "status": result["status"],
        "strict_DAE_l2": dae.get("native_final_residual_l2"),
        "state_DOFS": int(states.shape[1]),
        "steps": system.last_ts_steps,
        "heat_faces": heat_face,
        "wave_faces": wave_face,
        "physical_gates": gates,
        "max_total_CN_balance":
            physics["max_CN_step_H_balance_defect"],
        "max_interface_exchange_power":
            physics["max_individual_CN_interface_power"],
        "max_CN_interface_exchange_defect":
            physics["max_CN_interface_power_defect"],
    }, indent=2), flush=True)
    if not passed:
        raise SystemExit("Genuine CN completed, but energy physics gates failed")


if __name__ == "__main__":
    main()
