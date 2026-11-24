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
        report = system.audit_interface_ownership(10, (1, 2))
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
        "status": "PASS_DECLARED_DPHS_INTERFACE_OWNERSHIP_ONLY",
        "cases": reports,
        "original_example_source_unchanged": True,
        "getfem_weak_forms_reoriented": False,
        "PETSc_TS_executed": False,
        "initial_value_projection_executed": False,
        "physical_energy_qualified": False,
        "MPI2_qualified": False,
    }, indent=2))


if __name__ == "__main__":
    main()
