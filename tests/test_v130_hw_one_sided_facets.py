# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only guard on isolated one-sided 2D facet interpretation."""

import ast
from pathlib import Path
import unittest

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PROBE=ROOT/"sandbox/probe_v130_hw_one_sided_facets.py"
RUNNER=ROOT/"scripts/run-v130-hw-one-sided-facets.sh"


class OneSidedFacetContracts(unittest.TestCase):
    def test_one_sided_faces_are_from_2d_parent_region_only(self):
        s=PROBE.read_text(encoding="utf-8")
        for token in (
            "mesh.region(region)",
            "mesh.pid_from_cvid(candidates)",
            "candidates[node_counts==3]",
            "selection_based_on_geometric_connectivity_not_face_id",
            "bulk_region_geometric_classification",
            "mesh.outer_faces(2,volume)",
            "mesh.pid_in_faces(face)",
            "mesh.normal_of_faces(face)",
            "mesh.set_region(FACE_HEAT, disk)",
            "mesh.set_region(FACE_WAVE, inner)",
            "original_1d_region",
            "parent_face_normals",
        ):
            self.assertIn(token,s)

    def test_bulk_selector_accepts_mixed_face_flags_and_1d_convexes(self):
        """Execute the real isolated selector without GetFEM/PETSc.

        The old erroneous face==-1 assumption rejected this valid
        GetFEM region composition, although each triangle is a
        legitimate 2D convex parent.
        """
        tree=ast.parse(PROBE.read_text(encoding="utf-8"))
        function=next(
            node for node in tree.body
            if isinstance(node,ast.FunctionDef)
            and node.name=="bulk_2d_convexes"
        )
        compiled=compile(
            ast.fix_missing_locations(ast.Module(
                body=[function],type_ignores=[]
            )),
            str(PROBE),"exec",
        )
        scope={"np":np}
        exec(compiled,scope)

        class MixedRegion:
            def region(self,number):
                if number!=1:
                    raise ValueError(number)
                return np.array([
                    [10, 11, 12, 13],
                    [0, -1, 2, -1],
                ],dtype=np.int64)

            def pid_from_cvid(self,cvids):
                cells={
                    10:(1,2,3),   # 2D triangle, face tag 0
                    11:(2,3,4),   # 2D triangle, whole convex
                    12:(4,5),     # independent 1D element
                    13:(6,7),     # independent 1D element
                }
                pids=[]
                offsets=[0]
                for cv in cvids:
                    pids.extend(cells[int(cv)])
                    offsets.append(len(pids))
                return np.array(pids),np.array(offsets)

        volume, details=scope["bulk_2d_convexes"](MixedRegion(),1)
        self.assertEqual(volume.tolist(),[10,11])
        self.assertEqual(details["geometric_2D_triangle_parents"],2)
        self.assertEqual(details["excluded_1D_segment_convexes"],2)
        self.assertEqual(details["raw_face_entries"],2)
        self.assertEqual(details["raw_whole_convex_entries"],2)
        self.assertTrue(
            details["selection_based_on_geometric_connectivity_not_face_id"]
        )

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
        self.assertNotIn('Bulk volume must contain complete convexes',s)
        self.assertNotIn('np.any(ids[1,:]>=0)',s)

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
