# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Inspect original Heat-Wave DPHS declarations before initialization.

This reversible native diagnostic captures the fully declared ports,
bricks and controls at the FIRST initial-value call. It must not run
L2 projection, PETSc COO allocation or a time-integration step.
No original example is edited; no GetFEM weak form is reoriented.
"""

import json
import scrimp as S
from scrimp.dphs import DPHS
from examples.heat_wave import heat_wave_eq


class StopAtAudit(Exception):
    pass


def main():
    original_initial = DPHS.set_initial_value
    reports = []
    active = {"heat": None, "wave": None}

    def audit_before_initialization(system, *args, **kwargs):
        report = system.audit_interface_brick_plan(10, (1, 2))
        expected = {
            "Interface Heat": active["heat"],
            "Interface Wave": active["wave"],
        }
        if report["port_owners"] != expected:
            raise AssertionError(
                f"Original model's interface owners changed: {report['port_owners']}"
            )
        if report["facets_per_side"] != 78:
            raise AssertionError("Unexpected original Concentric interface")
        expected_bricks = {
            "M_Y_T": active["heat"],
            "B_T": active["heat"],
            "B_T^T": active["heat"],
            "M_Y_w": active["wave"],
            "B_w": active["wave"],
            "B_w^T": active["wave"],
            "M_U_T": active["heat"],
            "U_T_feedback": active["heat"],
            "M_U_w": active["wave"],
            "U_w_feedback": active["wave"],
        }
        actual_bricks = {
            entry["brick"]: entry["owner"]
            for entry in report["brick_plan"]
        }
        if actual_bricks != expected_bricks:
            raise AssertionError(
                f"Unexpected declared interface brick placement: {actual_bricks}"
            )
        if report["getfem_assembly_modified"]:
            raise AssertionError("Audit must not change actual assembly")
        reports.append({"heat_region": active["heat"],
                        "wave_region": active["wave"], **report})
        raise StopAtAudit

    try:
        DPHS.set_initial_value = audit_before_initialization
        for heat, wave in ((1, 2), (2, 1)):
            active.update(heat=heat, wave=wave)
            try:
                heat_wave_eq(
                    heat_region=heat, wave_region=wave,
                    final_time=.004, time_step=.001, save_step=.001,
                )
            except StopAtAudit:
                pass
    finally:
        DPHS.set_initial_value = original_initial

    if len(reports) != 2:
        raise AssertionError("Both original model permutations must be audited before projection")
    print(json.dumps({
        "schema": "scrimp-v130-native-declared-interface-ownership-v1",
        "status": "PASS_DECLARED_DPHS_INTERFACE_OWNERSHIP_AND_BRICK_PLAN_ONLY",
        "cases": reports,
        "original_example_source_unchanged": True,
        "getfem_weak_forms_reoriented": False,
        "all_ten_original_interface_bricks_placed_per_case": True,
        "PETSc_TS_executed": False,
        "initial_value_projection_executed": False,
        "physical_energy_qualified": False,
        "MPI2_qualified": False,
    }, indent=2))


if __name__ == "__main__":
    main()
