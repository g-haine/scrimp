#!/usr/bin/env python3
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Inspect real GetFEM shallow-water h/p VTK + original energy and GIF data.

No implicit assumption of exactly conservative BDF2, perfect spatial
accuracy or rotating viscous research validity. Only claim completion
and provenance of the original nonviscous FE visualization case.
"""
import argparse
from io import BytesIO
import json
from pathlib import Path
import sys
from xml.etree import ElementTree as ET

import numpy as np
from PIL import Image

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from examples.shallow_water import (  # noqa: E402
    true_pvd_files,native_triangular_field
)


def check(root):
    root=Path(root)
    result={
        "schema":"scrimp-v130-original-SWE-native-showcase-verify-v1",
        "status":"FAIL_OR_UNQUALIFIED_NATIVE_SWE_VISUALS",
        "errors":[],
        "original_FEM_and_PDE_unchanged":True,
        "BDF2_exact_H_conservation_qualified":False,
        "rotating_viscous_research_qualified":False,
        "MPI2_showcase_qualified":False,
        "spatial_convergence_qualified":False,
    }
    try:
        manifest=json.loads((root/"manifest.json").read_text())
        if (manifest.get("schema")!="v130-original-SWE-showcase-v1"
                or manifest.get("model")!="examples/shallow_water.py"
                or manifest.get("MPI")!=1
                or manifest.get("initialization")!="nonlinear_affine"
                or manifest.get("image_revision")!=
                   "1d6d0c90115213298047b28237a7f38087563a35"):
            raise RuntimeError("Wrong canonical model/image/native showcase manifest")
        native_rc=int((root/"native.exitcode").read_text())
        if native_rc:
            raise RuntimeError("Original GetFEM/PETSc BDF2 process failed")
        folder=root/"native"
        provenance=json.loads((folder/"provenance.json").read_text())
        if not (
            provenance.get("schema")==
                "v130-original-inviscid-shallow-water-native-visuals-v1"
            and provenance.get("source_of_spatial_values")==
                "original GetFEM h/p native PVD and VTU"
            and provenance.get("PDE_modified") is False
            and provenance.get("viscous_rotating_research_included") is False
            and provenance.get("H_exact_conservation_claimed") is False
            and provenance.get("initialization")==
                "nonlinear_affine (default direct example execution)"
            and provenance.get("velocity_definition")=="v=p/rho"
        ):
            raise RuntimeError("Not true canonical original GetFEM nonviscous fields")
        with np.load(folder/"physical_FE_traces.npz",allow_pickle=False) as npz:
            a={key:npz[key].copy() for key in npz.files}
        keys={"t","H","mass","kinetic","potential","boundary_power",
              "boundary_work","numerical_BDF2_energy_defect"}
        if set(a)!=keys:
            raise RuntimeError("Missing genuine original GetFEM quadratures")
        t=a["t"];H=a["H"]
        tf=float(manifest["final_time"])
        save=float(manifest["save_step"])
        if not (
            t.ndim==1 and len(t)>=6
            and all(x.shape==t.shape and np.all(np.isfinite(x))
                    for x in a.values())
            and abs(t[0])<1e-11 and abs(t[-1]-tf)<1e-8
            and np.all(np.diff(t)>0)
            and np.max(np.diff(t))<=save*1.001+1e-9
            and np.min(H)>0 and np.min(a["mass"])>0
            and np.allclose(H,a["kinetic"]+a["potential"],
                            atol=1e-6,rtol=1e-11)
            and np.allclose(a["numerical_BDF2_energy_defect"],
                            H-H[0]-a["boundary_work"],
                            atol=1e-6,rtol=1e-11)
            and abs(provenance["Hamiltonian_final"]-H[-1])<1e-8
            and provenance["true_FE_states"]==len(t)
        ):
            raise RuntimeError("Original FE time/H/mass or BDF2 defect inconsistent")
        if (np.max(abs(a["mass"]-a["mass"][0]))/a["mass"][0]>1e-5
                or np.max(abs(a["boundary_power"]))>1e-7):
            raise RuntimeError("Native mass or zero-control port work sanity failed")
        for field in ("h","p"):
            filenames=true_pvd_files(folder,field,t)
            for k in (0,len(t)//2,len(t)-1):
                xy,triangle,values=native_triangular_field(
                    filenames[k],field
                )
                if len(triangle)<10:
                    raise RuntimeError("Unexpected empty native rectangle")
                if field=="h" and not np.min(values)>0.:
                    raise RuntimeError("Negative original free surface")
                if field=="p" and not np.max(np.linalg.norm(values,axis=1))>=0.:
                    raise RuntimeError("Invalid native physical momentum")
            result[field+"_native_original_vtu"]=len(filenames)
        figs=folder/"figures"
        expected=[
            "saint_venant_H_masse_bilan.png",
            "saint_venant_H_masse_bilan.svg",
            "saint_venant_surface_vitesse.gif",
        ]
        for name in expected:
            file=figs/name
            if not file.is_file() or file.stat().st_size<300:
                raise RuntimeError("Missing in-example figure "+name)
            b=file.read_bytes()
            if name.endswith(".png") and b[:8]!=b"\x89PNG\r\n\x1a\n":
                raise RuntimeError("Bad in-example genuine PNG")
            if name.endswith(".gif") and b[:6] not in (b"GIF89a",b"GIF87a"):
                raise RuntimeError("Bad genuine native-state GIF")
            if name.endswith(".svg"):
                ET.fromstring(b)
        snapshots=sorted(figs.glob("saint_venant_champs_t*.png"))
        if len(snapshots)!=4:
            raise RuntimeError("Expected four native height/velocity panels")
        for f in snapshots:
            if f.read_bytes()[:8]!=b"\x89PNG\r\n\x1a\n":
                raise RuntimeError("Bad native height field snapshot")
        with Image.open(figs/"saint_venant_surface_vitesse.gif") as gif:
            framecount=gif.n_frames
            if not (5<=framecount<=36
                    and framecount==provenance["gif_original_FE_frames"]):
                raise RuntimeError("GIF frame count not genuine FE snapshots")
        selected=np.unique(np.linspace(
            0,len(t)-1,min(len(t),36),dtype=int
        ))
        if not np.allclose(
            np.asarray(provenance["original_selected_times"]),t[selected],
            atol=1e-8,rtol=0,
        ):
            raise RuntimeError("The GIF's recorded times are not original FE times")
        result.update({
            "status":"PASS_NATIVE_INVISCID_SWE_VISUALS",
            "source":manifest["source"],
            "MPI":1,
            "FE_states":int(len(t)),
            "original_H_initial":float(H[0]),
            "original_H_final":float(H[-1]),
            "max_relative_mass_drift":float(
                np.max(abs(a["mass"]-a["mass"][0]))/a["mass"][0]
            ),
            "final_signed_BDF2_energy_defect":float(
                a["numerical_BDF2_energy_defect"][-1]
            ),
            "GIF_frames_true_FE":framecount,
            "fixed_height_range":provenance["fixed_height_color_range"],
            "fixed_velocity_range":provenance["fixed_velocity_color_range"],
            "original_pvd_h_p_complete":True,
        })
    except Exception as exc:
        result["errors"].append(f"{type(exc).__name__}: {exc}")
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence",type=Path)
    args=parser.parse_args()
    root=args.evidence
    outcome=check(root)
    (root/"summary.json").write_text(
        json.dumps(outcome,indent=2)+"\n",encoding="utf-8"
    )
    print(json.dumps(outcome,indent=2),flush=True)
    raise SystemExit(0 if outcome["status"].startswith("PASS_") else 1)


if __name__=="__main__":
    main()
