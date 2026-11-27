# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Swapped heat=2/wave=1: compare historical and outward-normal laws.

Real GetFEM/PETSc stiffness assembly from the current Heat-Wave example.
The outward variant uses the source AS-IS. The historical variant
recreates the previous disk-based sign corrections *only inside this
sandbox diagnostic*, without editing the example.

A matrix contrast cannot by itself validate physical correctness.
No DAE/CN or energy claims.
"""

import argparse
import json
import os
from pathlib import Path

import numpy as np
from petsc4py import PETSc
from scipy.sparse import csr_matrix, save_npz

import scrimp as S
from scrimp.dphs import DPHS
from examples.heat_wave import heat_wave_eq


class EndAfterAssembly(Exception):
    """The native stiffness matrix was captured before DAE initialization."""


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--variant", choices=("published", "outward"), required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    if os.environ.get("SCRIMP_V130_AUTO_INTERFACE_PILOT") != "1":
        raise RuntimeError("Requires experimental automatic interface compiler")
    if PETSc.COMM_WORLD.getSize() != 1:
        raise RuntimeError("Swapped Heat-Wave sign diagnostic is MPI1 only")
    if args.out.exists() or args.out.with_suffix(".npz").exists():
        raise RuntimeError("Do not overwrite evidence files")
    if not args.out.parent.is_dir():
        raise RuntimeError("Parent output directory is not available")

    original_brick = S.Brick
    original_set_control = DPHS.set_control
    original_init = DPHS.init_step
    intercepted = {"brick": 0, "feedback": 0, "assembly": 0}

    def brick_constructor(*pos, **kwargs):
        if pos and pos[0] == "B_T":
            if len(pos) < 3 or pos[1] != "-U_T*Test_e_Q.Normal":
                raise RuntimeError("Current outward heat interface form changed")
            intercepted["brick"] += 1
            if args.variant == "published":
                pos = (pos[0], "U_T*Test_e_Q.Normal", *pos[2:])
        return original_brick(*pos, **kwargs)

    def control_law(system, name, expression):
        if name == "Interface Wave":
            if expression != "Y_T":
                raise RuntimeError("Current outward wave transmission law changed")
            intercepted["feedback"] += 1
            if args.variant == "published":
                expression = "-Y_T"
        return original_set_control(system, name, expression)

    def after_native_assembly(system):
        intercepted["assembly"] += 1
        if intercepted["assembly"] != 1:
            raise RuntimeError("Unexpected repeated actual assembly")
        compiler = system._internal_interface_compiler
        if compiler is None or not compiler.completed:
            raise RuntimeError("Original swapped model not auto-compiled")
        manifest = compiler.report["interfaces"]
        if len(manifest) != 1:
            raise RuntimeError("Expected precisely one logical interface")
        interface = manifest[0]
        if (interface["logical_region"] != 10
                or interface["port_owners"] != {
                    "Interface Heat": 2, "Interface Wave": 1
                } or interface["facet_count"] != 78):
            raise RuntimeError("Incorrect swapped physical subdomain owners")
        heat_port = system.ports["Interface Heat"]
        wave_port = system.ports["Interface Wave"]
        if (heat_port.get_region() != 10 or wave_port.get_region() != 10
                or heat_port.get_assembly_region() !=
                   interface["side_regions"][2]
                or wave_port.get_assembly_region() !=
                   interface["side_regions"][1]):
            raise RuntimeError("Incorrect distinct one-sided effective regions")

        expected_bt = ("-U_T*Test_e_Q.Normal"
                       if args.variant == "outward"
                       else "U_T*Test_e_Q.Normal")
        expected_control = ("(Y_T)*Test_U_w"
                            if args.variant == "outward"
                            else "(-Y_T)*Test_U_w")
        if (system.bricks["B_T"].get_form() != expected_bt
                or system.bricks["U_w_feedback"].get_form() != expected_control):
            raise RuntimeError("Original transmission sign variant not in model")

        n = int(system.gf_model.nbdof())
        offsets = {
            name: list(map(int, system.gf_model.interval_of_variable(name)))
            for name in ("T", "e_Q", "p", "q", "U_T", "Y_T",
                         "U_w", "Y_w", "U_bnd", "Y_bnd")
        }
        ptr, idx, data = system.stiffness.getValuesCSR()
        K = csr_matrix(
            (np.asarray(data).copy(), np.asarray(idx).copy(),
             np.asarray(ptr).copy()),
            shape=(n, n),
        )
        K.sum_duplicates()
        K.sort_indices()
        save_npz(args.out.with_suffix(".npz"), K, compressed=True)
        output = {
            "schema": "scrimp-v130-swapped-auto-interface-sign-K-v1",
            "status": "DIAGNOSTIC_SWAPPED_NATIVE_ASSEMBLY_ONLY",
            "variant": args.variant,
            "physical_meaning": (
                "pre-5042705 historical disk-based sign compensation, "
                "reconstructed only in the sandbox"
                if args.variant == "published"
                else "current outward-normal T=U_T and U_w=Y_T convention"
            ),
            "original_example_source_modified": False,
            "runtime_only_physical_sign_overrides": (
                ["B_T", "set_control(Interface Wave)"]
                if args.variant == "published" else []
            ),
            "logical_interface_region": 10,
            "effective_heat_region": heat_port.get_assembly_region(),
            "effective_wave_region": wave_port.get_assembly_region(),
            "compiler_manifest": compiler.report,
            "full_stiffness_shape": list(K.shape),
            "stiffness_stored_nnz": int(K.nnz),
            "intervals": offsets,
            "weak_heat_boundary_form": expected_bt,
            "wave_feedback_form": expected_control,
            "DAE_executed": False,
            "CN_energy_qualified": False,
            "MPI2_qualified": False,
        }
        args.out.write_text(json.dumps(output, indent=2) + "\n")
        print(json.dumps({
            "status": output["status"],
            "variant": args.variant,
            "shape": list(K.shape),
            "nnz": int(K.nnz),
            "heat_side": output["effective_heat_region"],
            "wave_side": output["effective_wave_region"],
        }, indent=2), flush=True)
        raise EndAfterAssembly

    S.Brick = brick_constructor
    DPHS.set_control = control_law
    DPHS.init_step = after_native_assembly
    try:
        try:
            heat_wave_eq(
                heat_region=2, wave_region=1,
                factor_solver="mumps",
                final_time=.004, time_step=.001, save_step=.001,
            )
        except EndAfterAssembly:
            pass
        else:
            raise RuntimeError("Actual native stiffness assembly was never intercepted")
    finally:
        S.Brick = original_brick
        DPHS.set_control = original_set_control
        DPHS.init_step = original_init

    if (intercepted != {"brick": 1, "feedback": 1, "assembly": 1}
            or not args.out.is_file()):
        raise RuntimeError(f"Unexpected model capture count: {intercepted}")


if __name__ == "__main__":
    main()
