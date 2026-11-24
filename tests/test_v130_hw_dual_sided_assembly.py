# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-only guards for reversible original Heat-Wave interface pilot."""

import ast
import importlib.util
from pathlib import Path
import unittest


ROOT=Path(__file__).resolve().parents[1]
PILOT=ROOT/"sandbox/probe_v130_hw_dual_sided_assembly.py"
RUN=ROOT/"scripts/run-v130-hw-dual-sided-assembly.sh"


class OriginalHeatWaveDualFaceAssemblyTests(unittest.TestCase):
    def test_unmodified_source_constructor_interception_is_restored(self):
        s=PILOT.read_text(encoding="utf-8")
        for token in (
            "original_Domain=S.Domain",
            "original_Control_Port=S.Control_Port",
            "original_Brick=S.Brick",
            "original_init=DPHS.init_step",
            "S.Domain=domain_wrapper",
            "S.Control_Port=control_wrapper",
            "S.Brick=brick_wrapper",
            "DPHS.init_step=before_native_DAE",
            "S.Domain=original_Domain",
            "S.Control_Port=original_Control_Port",
            "S.Brick=original_Brick",
            "DPHS.init_step=original_init",
            "raise AssemblyOnlyComplete",
        ):
            self.assertIn(token,s)

    def test_only_original_named_interface_ports_and_bricks_retargeted(self):
        s=PILOT.read_text(encoding="utf-8")
        for token in (
            '("M_Y_T", "B_T", "B_T^T")',
            '("M_Y_w", "B_w", "B_w^T")',
            '"Interface Heat": 30',
            '"Interface Wave": 31',
            "cells(old)!=(cells(disk)|cells(wave))",
            "new_one_sided_regions(mesh)",
            "compare_block_norms(",
            '"physics_energy_power_qualified":False',
            '"DAE_initialized":False',
        ):
            self.assertIn(token,s)

    def test_geometry_helper_is_standalone_and_explicitly_mounted(self):
        """Catch unmounted transitive sandbox dependencies before native run."""
        source=PILOT.read_text(encoding="utf-8")
        script=RUN.read_text(encoding="utf-8")
        helper=ROOT/"sandbox/v130_hw_interface_geometry.py"
        self.assertTrue(helper.is_file())
        self.assertIn(
            "from v130_hw_interface_geometry import new_one_sided_regions",
            source,
        )
        self.assertNotIn(
            "from probe_v130_hw_one_sided_facets import",
            source,
        )
        self.assertIn(
            '--volume "$root/sandbox/v130_hw_interface_geometry.py:',
            script,
        )
        self.assertNotIn(
            '--volume "$root/sandbox/probe_v130_hw_one_sided_facets.py:',
            script,
        )
        root_ast=ast.parse(helper.read_text(encoding="utf-8"))
        modules={
            alias.name.split(".")[0]
            for n in root_ast.body
            if isinstance(n,ast.Import)
            for alias in n.names
        }
        from_modules={
            n.module.split(".")[0] for n in root_ast.body
            if isinstance(n,ast.ImportFrom) and n.module
        }
        self.assertEqual(modules,{"numpy"})
        self.assertFalse(from_modules)
        self.assertEqual(
            {n.name for n in root_ast.body
             if isinstance(n,ast.FunctionDef)},
            {"region_description","bulk_2d_convexes","new_one_sided_regions"},
        )
        spec=importlib.util.spec_from_file_location(
            "v130_hw_interface_geometry",helper
        )
        module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertTrue(callable(module.new_one_sided_regions))

    def test_real_geometry_helper_executes_on_mixed_2d_mesh(self):
        """Exercise every helper global and region operation before Docker.

        The prior import-only test missed missing FACE_HEAT, FACE_WAVE
        and REGION_ORIGINAL globals. The synthetic P1 geometry has
        twelve disk faces, twelve annulus inner faces, four annulus
        outer faces, plus an original doubled region 10.
        """
        import numpy as np

        helper=ROOT/"sandbox/v130_hw_interface_geometry.py"
        spec=importlib.util.spec_from_file_location(
            "v130_hw_interface_geometry_runtime",helper
        )
        module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertEqual(
            (module.FACE_HEAT,module.FACE_WAVE,module.REGION_ORIGINAL),
            (30,31,10),
        )

        class P1ConcentricMesh:
            def __init__(self):
                self.data={
                    1:np.vstack((
                        np.arange(1,13,dtype=np.int64),
                        np.full(12,65535,dtype=np.int64),
                    )),
                    2:np.vstack((
                        np.arange(101,117,dtype=np.int64),
                        np.full(16,65535,dtype=np.int64),
                    )),
                    10:np.vstack((
                        np.r_[np.arange(1,13),np.arange(101,113)],
                        np.full(24,2,dtype=np.int64),
                    )),
                    20:np.vstack((
                        np.arange(113,117,dtype=np.int64),
                        np.full(4,2,dtype=np.int64),
                    )),
                }

            def region(self,region):
                return self.data[region]

            def regions(self):
                return np.asarray(tuple(self.data),dtype=np.int64)

            def pid_from_cvid(self,ids):
                ids=np.asarray(ids,dtype=np.int64).reshape(-1)
                points=np.stack((ids*3,ids*3+1,ids*3+2),axis=1)
                return points.reshape(-1),3*np.arange(len(ids)+1)

            def outer_faces(self,dim,ids):
                if dim!=2:
                    raise AssertionError("Expected parent triangles")
                ids=np.asarray(ids,dtype=np.int64)
                return np.vstack((ids,np.full(len(ids),2,dtype=np.int64)))

            def pid_in_faces(self,faces):
                cv=int(faces[0,0])
                return np.asarray((cv*10,cv*10+1),dtype=np.int64)

            def pts(self,points):
                cv=int(points[0]//10)
                r=0.6 if cv<113 else 1.0
                return np.asarray(((r,r),(0.,0.)))

            def normal_of_faces(self,faces):
                signs=np.where(faces[0]<100,1.,-1.)
                return np.vstack((signs,np.zeros(len(signs))))

            def set_region(self,region,faces):
                if region in self.data:
                    raise AssertionError("Must create new region only")
                self.data[region]=np.asarray(faces).copy()

        mesh=P1ConcentricMesh()
        original10=mesh.region(10).copy()
        topology=module.new_one_sided_regions(mesh)
        self.assertEqual(topology["disk_interface_2d_face_region"][
            "convex_face_pairs"],12)
        self.assertEqual(topology["annulus_interface_2d_face_region"][
            "convex_face_pairs"],12)
        self.assertEqual(topology["outer_annulus_face_count"],4)
        self.assertEqual(topology["original_1d_region"][
            "convex_face_pairs"],24)
        np.testing.assert_array_equal(mesh.region(10),original10)
        self.assertEqual(mesh.region(30).shape,(2,12))
        self.assertEqual(mesh.region(31).shape,(2,12))
        self.assertTrue(topology["no_volume_region_modified"])

    def test_runner_compares_both_original_assemblies_but_no_ts(self):
        s=RUN.read_text(encoding="utf-8")
        for token in (
            "for variant in baseline dual_sided",
            "--variant \"$variant\"",
            "SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo",
            "check-v130-published-immutable.py",
            "nonzero_bidirectional_feedback_blocks",
            "DAE_initialization_executed",
            "energy_or_physical_coupling_qualified",
        ):
            self.assertIn(token,s)


if __name__=="__main__":
    unittest.main()
