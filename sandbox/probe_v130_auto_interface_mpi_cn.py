# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Genuine MPI1/MPI2 automatic Heat-Wave strict DAE, CN and pH power.

Run the CURRENT outward-normal Heat-Wave example on real MPI ranks.
Only visual plotting and non-qualified parallel ParaView export are
suppressed; actual GetFEM/PETSc solve and all physical GetFEM
quadratures run on ALL ranks, not only rank zero.

Heat=1/wave=2 uses the independently qualified 'clamp' zero exterior
wave velocity initial-data candidate; Heat=2/wave=1 uses the original
Gaussian without a correction. The first changes only the prescribed
differential initial data, NOT the physical weak forms or solver.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
from mpi4py import MPI
import numpy as np
from petsc4py import PETSc

from scrimp.dphs import DPHS
from examples.heat_wave import heat_wave_eq
from measure_v130_hw_dual_facet_physics import (
    collect_heat_wave_dual_facet_physics,
)
from probe_v130_heat_wave_boundary_compatibility import (
    wave_boundary_partition,
)


def _real(array, label):
    values = np.asarray(array)
    if np.iscomplexobj(values):
        if not np.allclose(values.imag, 0., rtol=0., atol=1e-12):
            raise RuntimeError(f"Unexpected imaginary {label}")
        values = values.real
    values = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(values)):
        raise RuntimeError(f"Nonfinite {label}")
    return values


def _digest(array):
    return hashlib.sha256(
        np.ascontiguousarray(array, dtype=np.float64).tobytes()
    ).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--heat", type=int, choices=(1, 2), required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--pv-diagnostic", action="store_true",
        help="Isolate REAL original ParaView VTU/PVD exports per MPI rank; "
             "does not qualify shared production output",
    )
    args = parser.parse_args()
    wave = 3 - args.heat
    comm = MPI.COMM_WORLD
    rank, size = comm.rank, comm.size
    if (size not in (1, 2)
            or PETSc.COMM_WORLD.getSize() != size
            or os.environ.get("SCRIMP_V130_AUTO_INTERFACE_PILOT") != "1"):
        raise RuntimeError("Explicit real MPI1/MPI2 auto-interface gate required")
    if not args.out.is_dir():
        raise RuntimeError("Shared evidence output directory must exist")

    old_init = DPHS.init_step
    old_plot = DPHS.plot_Hamiltonian
    old_export = DPHS.export_to_pv
    old_show = plt.show
    skipped = {"plot": 0, "export": 0}
    init_capture = []

    def compatible_init(system):
        if init_capture:
            raise RuntimeError("Unexpected second strict DAE initialization")
        boundary = {
            "strategy": (
                "clamp-original-projected-p-on-wave-outer-boundary"
                if args.heat == 1 else "original-L2-Gaussian-no-adjustment"
            ),
            "initial_data_adjustment_applied": args.heat == 1,
        }
        if args.heat == 1:
            p0 = _real(system.gf_model.variable("p"), "projected p0").reshape(-1)
            T0 = _real(system.gf_model.variable("T"), "projected T0").copy()
            q0 = _real(system.gf_model.variable("q"), "projected q0").copy()
            B, F, _, _ = wave_boundary_partition(
                system.ports["p"].get_fem(), field_size=len(p0),
            )
            if len(p0) != 1105 or len(B) != 126:
                raise RuntimeError(
                    "Unexpected original CG1 exterior wave DOF layout"
                )
            corrected = p0.copy()
            corrected[B] = 0.
            if (not np.array_equal(corrected[F], p0[F])
                    or np.any(corrected[B])):
                raise RuntimeError("Wave boundary clamp corrupted FE coefficients")
            system.gf_model.set_variable(
                "p", np.ascontiguousarray(corrected, dtype=np.float64),
            )
            boundary.update({
                "boundary_wave_DOFs": int(len(B)),
                "wave_filtered_DOFs": int(len(p0)),
                "original_boundary_l2": float(np.linalg.norm(p0[B])),
                "corrected_boundary_l2": 0.,
            })
            if (not np.array_equal(
                    _real(system.gf_model.variable("T"), "T0"), T0)
                    or not np.array_equal(
                        _real(system.gf_model.variable("q"), "q0"), q0)):
                raise RuntimeError("Compatibility clamp modified heat/stress")
        old_init(system)  # genuine collective algebraic PETSc DAE initializer
        report = system.initialization_report
        if (not isinstance(report, dict) or report.get("status") != "PASS"
                or not report.get("differential_initial_exact")
                or report.get("pseudo_ts_executed")
                or float(report["native_final_residual_l2"]) >
                    float(report["feasibility_tolerance"])):
            raise RuntimeError("Strict DAE initial compatibility gate failed")
        init_capture.append(boundary)
        if rank == 0:
            print("STRICT_MPI_DAE_PASSED", json.dumps({
                "heat": args.heat, "ranks": size,
                "true_algebraic_residual": report["native_final_residual_l2"],
                "feasibility_tolerance": report["feasibility_tolerance"],
                "candidate": boundary["strategy"],
            }), flush=True)

    def skip_plot(system, *args_, **kwargs):
        skipped["plot"] += 1

    def skip_export(system, name_variable, path=None, t="All"):
        if args.pv_diagnostic:
            # Keep GetFEM's actual production export_to_pv() and all
            # collectives. Give every rank its own output directory:
            # production shared-file writes currently need a distinct
            # qualification gate because all ranks use the same names.
            rank_path = args.out / f"rank{rank}"
            target = rank_path / "pv" / name_variable
            if rank != 0:
                # Production exporter creates its leaf directory only
                # on rank 0. It must already exist on other ranks when
                # using an isolated path per rank.
                target.mkdir(parents=True, exist_ok=True)
            comm.Barrier()
            old_export(system, name_variable, path=str(rank_path), t=t)
            comm.Barrier()
        skipped["export"] += 1

    DPHS.init_step = compatible_init
    DPHS.plot_Hamiltonian = skip_plot
    DPHS.export_to_pv = skip_export
    plt.show = lambda *a, **k: None
    try:
        # Real live original GetFEM/PETSc model; only unsupported visual
        # I/O methods are intercepted, *never* FEM, Brick, Control, or TS.
        system = heat_wave_eq(
            heat_region=args.heat, wave_region=wave, factor_solver="mumps",
            final_time=0.004, time_step=0.001, save_step=0.001,
        )
    finally:
        DPHS.init_step = old_init
        DPHS.plot_Hamiltonian = old_plot
        DPHS.export_to_pv = old_export
        plt.show = old_show

    if len(init_capture) != 1 or skipped != {"plot": 1, "export": 2}:
        raise RuntimeError(
            f"Unexpected physical example execution count {init_capture}, {skipped}"
        )

    compiler = system._internal_interface_compiler
    if compiler is None or not compiler.completed:
        raise RuntimeError("No actual automatic interface compilation")
    manifest = compiler.report
    if len(manifest["interfaces"]) != 1:
        raise RuntimeError("Expected one physical interface")
    side = manifest["interfaces"][0]
    expected_ports = {"Interface Heat": args.heat, "Interface Wave": wave}
    if (side["logical_region"] != 10
            or side["facet_count"] != 78
            or side["port_owners"] != expected_ports):
        raise RuntimeError("Incorrect automatic oriented port ownership")
    heat_port = system.ports["Interface Heat"]
    wave_port = system.ports["Interface Wave"]
    heat_faces = int(heat_port.get_assembly_region())
    wave_faces = int(wave_port.get_assembly_region())
    if (heat_faces != side["side_regions"][args.heat]
            or wave_faces != side["side_regions"][wave]
            or heat_faces == wave_faces
            or heat_port.get_region() != 10
            or wave_port.get_region() != 10):
        raise RuntimeError("Inconsistent physical vs logical interface sides")

    times = _real(system.solution["t"], "accepted times").reshape(-1)
    states = _real(system.solution["z"], "full FE trajectory")
    n_expected = 7689 if args.heat == 1 else 8056
    if (times.shape != (5,) or not np.allclose(
            times, np.arange(5, dtype=float)*0.001, atol=1e-12, rtol=0.)
            or states.shape != (5, n_expected)
            or system.last_ts_converged_reason <= 0
            or system.last_ts_steps != 4
            or system.last_ts_rejections != 0
            or system.last_ts_snes_failures != 0):
        raise RuntimeError(
            f"MPI{size} physical 4-step PETSc CN trajectory invalid: "
            f"time={times.shape}, state={states.shape}, "
            f"steps={system.last_ts_steps}"
        )
    if rank == 0:
        print("MPI_AUTO_CN_COMPLETED", json.dumps({
            "heat": args.heat, "ranks": size, "steps": 4,
            "full_state_shape": list(states.shape),
        }), flush=True)

    # The original example has ALREADY evaluated get_Hamiltonian()
    # when preparing its logarithmic energy plot. Hamiltonian.compute()
    # appends values (it does not replace them), so invoking it a
    # second time would produce 10 entries for five accepted states.
    # Avoid that unintended recomputation on EVERY MPI rank; if the
    # current example ever changes to skip the calculation, perform
    # exactly one collective computation instead.
    already_computed = bool(system.hamiltonian.get_is_computed())
    if not already_computed:
        system.compute_Hamiltonian()
    term_lengths = {
        term.get_description(): len(term.get_values())
        for term in system.hamiltonian.get_terms()
    }
    if not term_lengths or any(
        count != len(times) for count in term_lengths.values()
    ):
        raise RuntimeError(
            f"Hamiltonian term history does not match the accepted "
            f"trajectory ({len(times)} entries): {term_lengths}"
        )

    # Port.compute remains a genuine GetFEM integration on ALL ranks.
    system.compute_powers()
    H = _real(system.get_Hamiltonian(), "physical Hamiltonian").reshape(-1)
    if H.shape != times.shape:
        raise RuntimeError("Hamiltonian time history has inconsistent length")
    traces, physical = collect_heat_wave_dual_facet_physics(
        system, times, H, heat_faces=heat_faces, wave_faces=wave_faces,
        heat_bulk=args.heat, wave_bulk=wave,
    )
    power_checks = {}
    for port, key in ((heat_port, "CN_P_interface_heat"),
                      (wave_port, "CN_P_interface_wave")):
        if port.get_power_sampling() != "cn_step":
            raise RuntimeError("Public Port.compute did not use CN step sampling")
        a = _real(port.get_power(), port.get_name()).reshape(-1)
        b = _real(traces[key], key).reshape(-1)
        if a.shape != b.shape or not np.allclose(a, b, atol=1e-10, rtol=1e-10):
            raise RuntimeError("Public interface power != independent sided trace")
        power_checks[port.get_name()] = float(np.max(np.abs(a-b)))

    scale = max(1., abs(float(H[0])))
    power_scale = max(1., float(physical["max_CN_heat_dissipation"]),
                      float(physical["max_individual_CN_interface_power"]))
    gates = {
        "strict_algebraic_DAE": True,
        "four_accepted_CN_steps": True,
        "nonzero_heat_dissipation":
            physical["max_CN_heat_dissipation"] > 1e-12,
        "nontrivial_interface_exchange":
            physical["interface_exchange_nontrivial"]
            and physical["max_individual_CN_interface_power"] > 1e-8,
        "total_CN_energy_balance":
            physical["max_CN_step_H_balance_defect"] <= 1e-7*scale,
        "heat_CN_energy_balance":
            physical["max_CN_step_heat_balance_defect"] <= 1e-7*scale,
        "wave_CN_energy_balance":
            physical["max_CN_step_wave_balance_defect"] <= 1e-7*scale,
        "closed_interface_exchange":
            physical["max_CN_interface_power_defect"] <= 1e-7*power_scale,
        "homogeneous_outer_boundary_work":
            physical["max_outer_port_power"] <= 1e-8*power_scale,
        "nonincreasing_H":
            physical["max_CN_positive_H_increment"] <= 1e-7*scale,
    }
    numerical = all(gates.values())

    # Every rank independently receives a full GetFEM model and a
    # complete TS trajectory in the current qualified code path.
    # Record fingerprints and check MPI-rank consensus independently
    # of comparing MPI1 and MPI2 runs offline.
    signature = {
        "times": _digest(times),
        "states": _digest(states),
        "H": _digest(H),
        "interface_heat_power": _digest(traces["CN_P_interface_heat"]),
        "interface_wave_power": _digest(traces["CN_P_interface_wave"]),
    }
    signatures = comm.allgather(signature)
    rank_consensus = all(item == signatures[0] for item in signatures)
    local_result = {
        "schema": "scrimp-v130-native-auto-interface-MPI-CN-v1",
        "status": (
            "PASS_RANK_MPI_CN_AND_PHYSICAL_GATES"
            if numerical and rank_consensus
            else "FAIL_OR_UNQUALIFIED_MPI_CN_PHYSICS"
        ),
        "heat_region": args.heat, "wave_region": wave,
        "MPI_ranks": size, "rank": rank,
        "initial_candidate": init_capture[0],
        "true_strict_DAE_initializer": system.initialization_report,
        "hamiltonian_was_already_computed_by_example": already_computed,
        "hamiltonian_terms_have_one_value_per_accepted_time": True,
        "hamiltonian_term_lengths": term_lengths,
        "real_PETSc_CN_steps": int(system.last_ts_steps),
        "full_state_shape": list(states.shape),
        "automatic_interface": manifest,
        "heat_side_region": heat_faces,
        "wave_side_region": wave_faces,
        "rank_consensus": rank_consensus,
        "rank_fingerprints": signatures,
        "physical_gates": gates,
        "physical_diagnostics": physical,
        "public_port_power_max_differences": power_checks,
        "plot_and_ParaView_export_skipped_on_all_ranks":
            not args.pv_diagnostic,
        "MPI_rank_isolated_real_ParaView_export_diagnostic":
            args.pv_diagnostic,
        "parallel_PNG_and_ParaView_shared_export_qualified": False,
        "physical_GetFEM_assembly_executed_on_all_ranks": True,
        "source_model_signs_unchanged": True,
        "MPI2_parallel_visual_exports_qualified": False,
    }
    np.savez_compressed(
        args.out/f"rank{rank}.npz",
        times=times, states=states, H=H,
        **{"phys_"+key: x for key, x in traces.items()},
    )
    (args.out/f"rank{rank}.json").write_text(
        json.dumps(local_result, indent=2)+"\n", encoding="utf-8"
    )
    comm.Barrier()
    if rank == 0:
        print("MPI_AUTO_CN_PHYSICS", json.dumps({
            "status": local_result["status"], "heat": args.heat,
            "ranks": size, "rank_consensus": rank_consensus,
            "max_total_energy_defect":
                physical["max_CN_step_H_balance_defect"],
            "max_interface_power_defect":
                physical["max_CN_interface_power_defect"],
        }), flush=True)
    if not numerical or not rank_consensus:
        raise RuntimeError("MPI CN/physical gate failed; preserved native evidence")


if __name__ == "__main__":
    main()
