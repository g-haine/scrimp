# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Experimental automatic Heat-Wave: real four-step strict-DAE/CN physics.

NO override to published Domain, Control_Port, Brick or heat_wave.py.
Uses SCRIMP_V130_AUTO_INTERFACE_PILOT=1 and the original published
example / DAE-compatible initial-data candidate. Read the effective
one-sided GetFEM regions from the compiled, original logical interface
ports, instead of assuming historical fixed region IDs 30 and 31.

This pilot deliberately does not claim MPI2, swapped-region physical
sign convention, or fully generic nonconforming interface support.
"""

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

import probe_v130_heat_wave_compatible_cn as compatible_cn
import probe_v130_main_2d_examples as main2d
import measure_v130_hw_dual_facet_physics as facet_physics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", choices=("clamp", "constrained_l2"),
                        required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    if MPI.COMM_WORLD.size != 1 or PETSc.COMM_WORLD.getSize() != 1:
        raise RuntimeError("Automatic Heat-Wave CN pilot is MPI1 only")
    if os.environ.get("SCRIMP_V130_AUTO_INTERFACE_PILOT") != "1":
        raise RuntimeError("Automatic interface pilot must be explicitly enabled")
    if args.out.exists() and any(args.out.iterdir()):
        raise RuntimeError("Refusing existing evidence directory")

    old_collector = main2d.collect_heat_wave_physics
    old_heat = facet_physics.HEAT_FACES
    old_wave = facet_physics.WAVE_FACES
    old_argv = sys.argv
    captures = []

    def collect_auto_physics(system, times, hamiltonian):
        if captures:
            raise RuntimeError("Duplicate actual physical collector")
        compiler = system._internal_interface_compiler
        if compiler is None or not compiler.completed or compiler.report is None:
            raise RuntimeError("Automatic DPHS interface compilation did not run")
        report = compiler.report
        if len(report["interfaces"]) != 1:
            raise RuntimeError("Expected exactly one automatically compiled interface")
        interface = report["interfaces"][0]
        if (interface["logical_region"] != 10
                or interface["subdomains"] != [1, 2]
                or interface["facet_count"] != 78
                or interface["port_owners"] != {
                    "Interface Heat": 1, "Interface Wave": 2
                }):
            raise RuntimeError("Original logical interface ownership changed")

        heat_port = system.ports["Interface Heat"]
        wave_port = system.ports["Interface Wave"]
        if heat_port.get_region() != 10 or wave_port.get_region() != 10:
            raise RuntimeError("User-visible logical interface was overwritten")
        heat = int(heat_port.get_assembly_region())
        wave = int(wave_port.get_assembly_region())
        if heat == wave or heat == 10 or wave == 10:
            raise RuntimeError("One-sided interface regions were not resolved")
        if heat != int(interface["side_regions"][1]):
            raise RuntimeError("Thermal effective region differs from compiler")
        if wave != int(interface["side_regions"][2]):
            raise RuntimeError("Wave effective region differs from compiler")

        original_mesh_region = np.asarray(
            system.domain.get_mesh()[0].region(10), dtype=np.int64
        )
        if original_mesh_region.shape != (2, 156):
            raise RuntimeError("Original doubled logical region was overwritten")

        # Existing physical balance collector is kept byte-for-byte
        # unchanged. Its two region identifiers are replaced ONLY for
        # the duration of this actual physical midpoint measurement.
        # All numerical definitions and thresholds remain identical.
        facet_physics.HEAT_FACES = heat
        facet_physics.WAVE_FACES = wave
        try:
            traces, physical = facet_physics.collect_heat_wave_dual_facet_physics(
                system, times, hamiltonian,
            )
        finally:
            facet_physics.HEAT_FACES = old_heat
            facet_physics.WAVE_FACES = old_wave
        if (physical["interface_heat_region"] != heat
                or physical["interface_wave_region"] != wave
                or not physical["old_doubled_interface_region_excluded_from_physical_power"]):
            raise RuntimeError("Physical power used incorrect interface faces")

        # Qualify the PUBLIC algebraic Port.compute API as well as the
        # separate GetFEM physical collector. The port must use its
        # resolved effective side; evaluating on doubled region 10
        # would generally give the wrong interface power.
        port_power_defects = {}
        for port, trace_name in (
            (heat_port, "CN_P_interface_heat"),
            (wave_port, "CN_P_interface_wave"),
        ):
            port.compute(
                system.solution, system.gf_model, system.domain,
                sampling="cn_step",
            )
            actual = np.asarray(port.get_power(), dtype=float).reshape(-1)
            expected = np.asarray(traces[trace_name], dtype=float).reshape(-1)
            if actual.shape != expected.shape or not np.all(np.isfinite(actual)):
                raise RuntimeError("Invalid publicly integrated interface port power")
            discrepancy = float(np.max(np.abs(actual - expected)))
            if not np.allclose(actual, expected, rtol=1e-10, atol=1e-10):
                raise RuntimeError(
                    f"Port.compute disagrees with sided physical trace "
                    f"for {port.get_name()}: {discrepancy}"
                )
            port_power_defects[port.get_name()] = discrepancy
        physical["public_port_compute_uses_resolved_side"] = True
        physical["public_port_compute_cn_max_defect"] = port_power_defects

        captures.append({
            "compiler": report,
            "one_sided_heat_region": heat,
            "one_sided_wave_region": wave,
            "original_logical_region": 10,
            "original_logical_face_entries": int(original_mesh_region.shape[1]),
            "one_sided_power_measurement": True,
        })
        return traces, physical

    main2d.collect_heat_wave_physics = collect_auto_physics
    try:
        sys.argv = [
            Path(__file__).name,
            "--candidate", args.candidate,
            "--backend", "mumps",
            "--out", str(args.out),
        ]
        compatible_cn.main()
    finally:
        main2d.collect_heat_wave_physics = old_collector
        facet_physics.HEAT_FACES = old_heat
        facet_physics.WAVE_FACES = old_wave
        sys.argv = old_argv

    if len(captures) != 1:
        raise RuntimeError("Original CN physical postprocessing was not reached")
    if not (args.out / "compatible_cn_report.json").is_file():
        raise RuntimeError("Original compatible strict-DAE/CN report missing")
    dae = json.loads((args.out / "compatible_cn_report.json").read_text())
    original = json.loads((args.out / "original_2d_example_report.json").read_text())
    physical = json.loads((args.out / "heat_wave_physics_diagnostic.json").read_text())
    if (not dae["true_strict_DAE_pre_Ts"]
            or dae["DAE_correction"]["DAE_initializer"]["status"] != "PASS"
            or not dae["DAE_correction"]["DAE_initializer"][
                "differential_initial_exact"
            ]
            or original["steps"] != 4
            or original["rejected_steps"] != 0
            or original["SNES_failures"] != 0
            or original["PETSc_reason"] <= 0
            or original["full_state_DOFs"] != 7689
            or not physical["interface_exchange_nontrivial"]):
        raise RuntimeError("Original strict DAE/CN or nontrivial power gates failed")

    science = bool(dae["physical_balance_gate_passed"])
    result = {
        "schema": "scrimp-v130-experimental-auto-interface-real-CN-v1",
        "status": (
            "PASS_AUTO_CN_BALANCE_AND_NONTRIVIAL_PORT_EXCHANGE"
            if science else
            "DIAGNOSTIC_AUTO_CN_CONVERGED_PHYSICAL_BALANCE_OPEN"
        ),
        "candidate": args.candidate,
        "source": "unchanged published heat_wave_eq",
        "no_constructor_region_remapping": True,
        "automatic_interface_manifest": captures[0],
        "actual_true_strict_DAE_initialization_passed": True,
        "actual_real_4_step_PETSc_CN_passed": True,
        "discrete_CN_energy_gate_passed": science,
        "nontrivial_interface_exchange_passed":
            bool(physical["interface_exchange_nontrivial"]),
        "actual_physical_diagnostics": physical,
        "MPI2_qualified": False,
        "region_swapped_physical_model_qualified": False,
        "automatic_solver_policy_changed": False,
    }
    target = args.out / "auto_interface_cn_audit.json"
    if target.exists():
        raise RuntimeError("Refuse to overwrite original auto CN report")
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "status": result["status"],
        "candidate": args.candidate,
        "steps": original["steps"],
        "one_sided_heat_region": captures[0]["one_sided_heat_region"],
        "one_sided_wave_region": captures[0]["one_sided_wave_region"],
        "DAE_residual": dae["DAE_correction"]["DAE_initializer"][
            "native_final_residual_l2"
        ],
        "max_CN_total_balance_defect":
            physical["max_CN_step_H_balance_defect"],
        "max_CN_interface_power_defect":
            physical["max_CN_interface_power_defect"],
        "max_individual_CN_interface_power":
            physical["max_individual_CN_interface_power"],
        "interface_excitation_relative_to_heat_dissipation":
            physical["interface_exchange_ratio_to_heat_dissipation"],
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
