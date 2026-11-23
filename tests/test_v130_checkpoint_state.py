# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.

"""#29 host-only checkpoint schema, mesh, physics, corruption, MPI guards.

No synthetic success in this file is native PETSc/GetFEM evidence.
"""

import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"sandbox"))

from v130_checkpoint_state import (
    canonical, collective_read, collective_write, model_provenance,
    read_checkpoint, sha256, write_checkpoint,
)


class Mesh:
    def __init__(self, changed=False):
        self.x = np.array([[0., 1., 0.], [0., 0., 1.]])
        if changed:
            self.x[1, 2] = 1.01

    def pts(self): return self.x
    def dim(self): return 2
    def cvid(self): return np.array([0], dtype=np.int32)
    def pts_from_cvid(self):
        return self.x.copy(), np.array([0,3], dtype=np.int32)


class Fem:
    def basic_dof_nodes(self):
        return np.array([[0., 1., 0.], [0., 0., 1.]])
    def nbdof(self): return 3


class Parameter:
    def __init__(self,value): self.value=value
    def get_name(self): return "g"
    def get_kind(self): return "scalar"
    def get_expression(self): return self.value


class Port:
    def __init__(self, p): self.p=p
    def get_name(self): return "h"
    def get_flow(self): return "h"
    def get_effort(self): return "e_h"
    def get_kind(self): return "scalar-field"
    def get_region(self): return None
    def get_mesh_id(self): return 0
    def get_algebraic(self): return False
    def get_substituted(self): return False
    def get_fem(self): return Fem()
    def get_parameters(self): return [Parameter(self.p)]


class Brick:
    def __init__(self,control): self.control=control
    def get_form(self): return f"h*Test_h + ({self.control})"
    def get_position(self): return "effort"
    def get_dt(self): return False
    def get_linear(self): return True
    def get_explicit(self): return False
    def get_regions(self): return [1]
    def get_mesh_id(self): return 0


class Domain:
    def __init__(self,mesh): self.mesh=mesh
    def get_mesh(self): return [self.mesh]
    def get_name(self): return "unit triangle"
    def get_dim(self): return [2]
    def get_subdomains(self): return [{"Omega":1}]
    def get_boundaries(self): return [{"Gamma":10}]


class Model:
    def nbdof(self): return 6
    def interval_of_variable(self,name):
        return {"h":(0,3), "e_h":(3,3)}[name]


class System:
    def __init__(self,changed_mesh=False,control="0.",p="0.1"):
        self.domain=Domain(Mesh(changed_mesh))
        self.gf_model=Model()
        self.ports={"h":Port(p)}
        self.bricks={"weak":Brick(control)}


class MPI1:
    rank=0
    def allgather(self,v): return [v]
    def bcast(self,v,root=0): return v


class MPI2Mismatch(MPI1):
    def allgather(self,v): return [v,("different-state", "other", .03, True)]


class CheckpointSchemaTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.file=Path(self.temp.name)/"checkpoint.npz"
        self.z=np.array([5.,5.,5.,.5,.5,.5])
        self.provenance=model_provenance(System())

    def save(self):
        return write_checkpoint(
            self.file,z=self.z,t=.02,
            provenance=self.provenance,ts_type="cn",
            accepted=True,extra={"origin":"synthetic-host-test"}
        )

    def test_real_atomic_roundtrip_and_provenance(self):
        m=self.save()
        t,z,header=read_checkpoint(
            self.file,current_provenance=self.provenance
        )
        self.assertEqual(t,.02)
        np.testing.assert_array_equal(z,self.z)
        self.assertEqual(header,m)
        self.assertEqual(list(Path(self.temp.name).glob("*.pending")),[])

    def test_mesh_connectivity_fingerprint_rejects_change(self):
        self.save()
        other=model_provenance(System(changed_mesh=True))
        with self.assertRaisesRegex(ValueError,"mesh, FE numbering"):
            read_checkpoint(self.file,current_provenance=other)

    def test_physical_boundary_law_must_match(self):
        self.save()
        other=model_provenance(System(control="omega*Tangent"))
        with self.assertRaisesRegex(ValueError,"mesh, FE numbering"):
            read_checkpoint(self.file,current_provenance=other)

    def test_physical_parameter_must_match(self):
        self.save()
        other=model_provenance(System(p=".101"))
        with self.assertRaisesRegex(ValueError,"mesh, FE numbering"):
            read_checkpoint(self.file,current_provenance=other)

    def test_no_unaccepted_source_state(self):
        with self.assertRaisesRegex(ValueError,"non-accepted"):
            write_checkpoint(self.file,z=self.z,t=.02,
                provenance=self.provenance,ts_type="cn",accepted=False)
        self.assertFalse(self.file.exists())

    def test_no_nonfinite_coefficients(self):
        z=self.z.copy();z[2]=float("nan")
        with self.assertRaisesRegex(ValueError,"NaN or infinite"):
            write_checkpoint(self.file,z=z,t=.02,
                provenance=self.provenance,ts_type="cn",accepted=True)

    def test_no_missing_algebraic_components(self):
        with self.assertRaisesRegex(ValueError,"DOF count mismatch"):
            write_checkpoint(self.file,z=self.z[:3],t=.02,
                provenance=self.provenance,ts_type="cn",accepted=True)

    def test_no_silent_bdf_history_restart(self):
        with self.assertRaisesRegex(ValueError,"only PETSc TSCN"):
            write_checkpoint(self.file,z=self.z,t=.02,
                provenance=self.provenance,ts_type="bdf",accepted=True)
        self.save()
        with self.assertRaisesRegex(ValueError,"requires original PETSc CN"):
            read_checkpoint(self.file,current_provenance=self.provenance,ts_type="bdf")

    def test_tampered_state_hash_rejected(self):
        self.save()
        with np.load(self.file,allow_pickle=False) as x:
            state=x["state"].copy()
            manifest=x["manifest"].copy()
        state[0]+=1.
        np.savez_compressed(self.file,state=state,manifest=manifest)
        with self.assertRaisesRegex(ValueError,"Damaged checkpoint FE state"):
            read_checkpoint(self.file,current_provenance=self.provenance)

    def test_tampered_metadata_model_hash_rejected(self):
        self.save()
        with np.load(self.file,allow_pickle=False) as x:
            state=x["state"].copy()
            header=json.loads(x["manifest"].tobytes().decode("utf-8"))
        header["model"]["detail"]["bricks"][0]["expression"]="different"
        np.savez_compressed(
            self.file,state=state,
            manifest=np.frombuffer(canonical(header),dtype=np.uint8)
        )
        with self.assertRaisesRegex(ValueError,"mesh, FE numbering"):
            read_checkpoint(self.file,current_provenance=self.provenance)

    def test_truncated_file_fails_closed(self):
        self.save()
        with self.file.open("rb+") as f:
            f.truncate(50)
        with self.assertRaises(Exception):
            read_checkpoint(self.file,current_provenance=self.provenance)

    def test_collective_single_rank_roundtrip(self):
        mpi=MPI1()
        collective_write(self.file,mpi=mpi,z=self.z,t=.02,
                         provenance=self.provenance,ts_type="cn",
                         accepted=True)
        t,z,_=collective_read(
            self.file,mpi=mpi,current_provenance=self.provenance
        )
        self.assertEqual(t,.02)
        np.testing.assert_array_equal(z,self.z)

    def test_collective_mpi_state_mismatch_rejected(self):
        with self.assertRaisesRegex(ValueError,"MPI inconsistent"):
            collective_write(
                self.file,mpi=MPI2Mismatch(),z=self.z,t=.02,
                provenance=self.provenance,ts_type="cn",accepted=True
            )

    def test_corrupted_schema_rejected(self):
        self.save()
        with np.load(self.file,allow_pickle=False) as x:
            state=x["state"].copy()
            header=json.loads(x["manifest"].tobytes().decode("utf-8"))
        header["schema"]="another"
        np.savez_compressed(
            self.file,state=state,
            manifest=np.frombuffer(canonical(header),dtype=np.uint8)
        )
        with self.assertRaisesRegex(ValueError,"Unsupported checkpoint schema"):
            read_checkpoint(self.file,current_provenance=self.provenance)


if __name__=="__main__":
    unittest.main()
