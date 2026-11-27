# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Original, unchanged Heat-Wave with experimental automatic interface compile.

Capture the full GetFEM/PETSc K matrix before DAE initialization or TS.
Source-only experimental mechanism SCRIMP_V130_AUTO_INTERFACE_PILOT=1.
No manual Domain, Control_Port, Brick or geometry monkey-patching.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix, save_npz
from petsc4py import PETSc

from scrimp.dphs import DPHS
from examples.heat_wave import heat_wave_eq


class AssemblyOnlyComplete(Exception):
    """Expected exit after the complete, actual native matrix assembly."""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if os.environ.get("SCRIMP_V130_AUTO_INTERFACE_PILOT") != "1":
        raise RuntimeError("Experimental automatic compiler must be enabled")
    if PETSc.COMM_WORLD.getSize() != 1:
        raise RuntimeError("Pilot is MPI1 only until validated")
    if args.out.exists() or args.out.with_suffix(".npz").exists():
        raise RuntimeError("Refuse to overwrite native evidence")
    if not args.out.parent.is_dir():
        raise RuntimeError("Output directory does not exist")

    original = DPHS.init_step
    captured = []

    def after_actual_assembly(system):
        if captured:
            raise RuntimeError("Duplicate assembly interception")
        compiler = system._internal_interface_compiler
        if compiler is None or not compiler.completed:
            raise RuntimeError("Internal interface was not automatically compiled")
        report = compiler.report
        if len(report["interfaces"]) != 1:
            raise RuntimeError("Expected exactly one discovered interface")
        interface = report["interfaces"][0]
        if interface["logical_region"] != 10:
            raise RuntimeError("Unexpected user logical interface")
        if interface["port_owners"] != {
            "Interface Heat": 1, "Interface Wave": 2
        }:
            raise RuntimeError("Unexpected automatic port ownership")
        if interface["facet_count"] != 78:
            raise RuntimeError("Unexpected number of interface facets")
        expected = {
            "M_Y_T": 1, "B_T": 1, "B_T^T": 1,
            "M_Y_w": 2, "B_w": 2, "B_w^T": 2,
            "M_U_T": 1, "U_T_feedback": 1,
            "M_U_w": 2, "U_w_feedback": 2,
        }
        if interface["brick_owners"] != expected:
            raise RuntimeError("Unexpected auto-selected brick owners")
        if any(
            system.ports[name].get_region() != 10
            for name in ("Interface Heat", "Interface Wave")
        ):
            raise RuntimeError("The public interface region was overwritten")
        mesh = system.domain.get_mesh()[0]
        original = np.asarray(mesh.region(10), dtype=np.int64)
        side_regions = interface["side_regions"]
        for bulk in (1, 2):
            region = int(side_regions[bulk])
            selected = np.asarray(mesh.region(region), dtype=np.int64)
            expected_side = system.domain.resolve_interface(10, (1, 2)).for_subdomain(bulk)
            if set(map(tuple, selected.T)) != set(map(tuple, expected_side.T)):
                raise RuntimeError("Compiled facet region differs from topology")
        n = int(system.gf_model.nbdof())
        if n != 7689:
            raise RuntimeError(f"Unexpected real GetFEM DOF count: {n}")
        ptr, idx, values = system.stiffness.getValuesCSR()
        K = csr_matrix(
            (np.asarray(values).copy(), np.asarray(idx).copy(),
             np.asarray(ptr).copy()),
            shape=(n, n),
        )
        K.sum_duplicates()
        K.sort_indices()
        save_npz(args.out.with_suffix(".npz"), K, compressed=True)
        names = ("T", "p", "q", "e_Q", "U_T", "Y_T",
                 "U_w", "Y_w", "U_bnd", "Y_bnd")
        intervals = {
            name: list(map(int, system.gf_model.interval_of_variable(name)))
            for name in names
        }
        output = {
            "schema": "scrimp-v130-experimental-production-path-auto-K-v1",
            "status": "DIAGNOSTIC_ACTUAL_AUTO_GETFEM_K_ONLY",
            "compiler": report,
            "full_stiffness_shape": list(K.shape),
            "stiffness_nnz": int(K.nnz),
            "variable_intervals": intervals,
            "original_user_interface_region": 10,
            "original_interface_face_entries": int(original.shape[1]),
            "original_example_source_unchanged": True,
            "no_constructor_remapping": True,
            "getfem_assembly_executed": True,
            "DAE_initialized": False,
            "CN_or_energy_qualified": False,
            "MPI2_qualified": False,
        }
        args.out.write_text(json.dumps(output, indent=2) + "\n")
        captured.append(output)
        print("EXPERIMENTAL_FULL_NATIVE_AUTO_K_ASSEMBLED", K.shape, K.nnz,
              flush=True)
        raise AssemblyOnlyComplete

    DPHS.init_step = after_actual_assembly
    try:
        try:
            heat_wave_eq(
                heat_region=1, wave_region=2, factor_solver="mumps",
                final_time=.004, time_step=.001, save_step=.001,
            )
        except AssemblyOnlyComplete:
            pass
        else:
            raise RuntimeError("The native assembly interception was not reached")
    finally:
        DPHS.init_step = original
    if len(captured) != 1 or not args.out.is_file():
        raise RuntimeError("Actual auto K was not captured exactly once")


if __name__ == "__main__":
    main()
