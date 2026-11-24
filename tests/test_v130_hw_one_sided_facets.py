# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only guard on isolated one-sided 2D facet interpretation."""

from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
PROBE=ROOT/"sandbox/probe_v130_hw_one_sided_facets.py"
RUNNER=ROOT/"scripts/run-v130-hw-one-sided-facets.sh"


class OneSidedFacetContracts(unittest.TestCase):
    def test_one_sided_faces_are_from_2d_parent_region_only(self):
        s=PROBE.read_text(encoding="utf-8")
        for token in (
            "mesh.region(number)",
            "mesh.outer_faces(2,volume)",
            "mesh.pid_in_faces(face)",
            "mesh.normal_of_faces(face)",
            "mesh.set_region(FACE_HEAT, disk)",
            "mesh.set_region(FACE_WAVE, inner)",
            "original_1d_region",
            "parent_face_normals",
        ):
            self.assertIn(token,s)

    def test_no_model_mutation_and_native_energy_collector_first(self):
        s=PROBE.read_text(encoding="utf-8")
        for token in (
            "traces,physics=original_collector(system,times,H,*a,**kw)",
            "system.get_quantity(",
            "original.collect_heat_wave_physics=observer",
            "original.collect_heat_wave_physics=original_collector",
            "original_source_and_weak_forms_unmodified",
            "physical_interface_coupling_qualified",
            "unavailable_cross_FEM_traces",
        ):
            self.assertIn(token,s)
        self.assertNotIn("system.stiffness.zeroEntries",s)

    def test_native_runner_is_pinned_and_readonly(self):
        s=RUNNER.read_text(encoding="utf-8")
        for token in (
            "check-v130-published-immutable.py",
            "SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo",
            "probe_v130_hw_one_sided_facets.py",
            "PASS_FACET_DIAGNOSTIC_CAPTURED_PHYSICS_REQUIRES_FOLLOWUP",
            "mpiexec --oversubscribe --bind-to none -n 1",
        ):
            self.assertIn(token,s)


if __name__=="__main__":
    unittest.main()
