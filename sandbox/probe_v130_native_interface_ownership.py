# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Inspect original Heat-Wave DPHS declarations before any TS time step.

This is a reversible native diagnostic of the current production API:
no source example is edited; no GetFEM matrix is reassembled or
reoriented by the new ownership analysis. Exit deliberately at solve().
"""

import json
import scrimp as S
from scrimp.dphs import DPHS
from examples.heat_wave import heat_wave_eq


class StopAtAudit(Exception):
    pass


def main():
    original_solve = DPHS.solve
    reports = []

    def audit_at_solve(system):
        report = system.audit_interface_ownership(10, (1, 2))
        expected = {"Interface Heat": 1, "Interface Wave": 2}
        if report["port_owners"] != expected:
            raise AssertionError(
                f"Original model's interface owners changed: {report['port_owners']}"
            )
        if report["facets_per_side"] != 78:
            raise AssertionError("Unexpected original Concentric interface")
        reports.append(report)
        raise StopAtAudit

    try:
        DPHS.solve = audit_at_solve
        try:
            heat_wave_eq(
                heat_region=1, wave_region=2,
                final_time=.004, time_step=.001, save_step=.001,
            )
        except StopAtAudit:
            pass
    finally:
        DPHS.solve = original_solve

    if len(reports) != 1:
        raise AssertionError("Original model did not reach ownership audit")
    print(json.dumps({
        "schema": "scrimp-v130-native-declared-interface-ownership-v1",
        "status": "PASS_DECLARED_DPHS_INTERFACE_OWNERSHIP_ONLY",
        "report": reports[0],
        "original_example_source_unchanged": True,
        "getfem_weak_forms_reoriented": False,
        "PETSc_TS_executed": False,
        "physical_energy_qualified": False,
        "MPI2_qualified": False,
    }, indent=2))


if __name__ == "__main__":
    main()
