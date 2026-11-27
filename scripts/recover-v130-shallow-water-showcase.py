#!/usr/bin/env python3
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Rescue REAL CG2 h / DG1 p spatial figures from a completed native VTK ZIP.

The original source-14181a94 1200-step MUMPS BDF2 simulation converged
at t=0.12 and exported all 61 h/p VTK snapshots before an OLD renderer
rejected valid type-22 quadratic triangles. The VTK files are original,
not synthetic, and can be rendered without repeating the 1200-step solve.

This command NEVER reclassifies the original ZIP's failed native renderer.
Original H/mass quadratures are NOT in that ZIP: they cannot be inferred
from the saved Float32 nodes; no energy plot or physical gate is claimed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys
from tempfile import TemporaryDirectory
from xml.etree import ElementTree as ET
from zipfile import ZipFile, ZIP_DEFLATED

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from examples.shallow_water import (  # noqa: E402
    render_original_shallow_water_spatial_VTU,
)


def recover(archive, output, *, frames=36):
    archive=Path(archive).resolve()
    output=Path(output).resolve()
    if output.exists():
        if any(output.iterdir()):
            raise FileExistsError("Refusing overwrite of existing native recovery")
    else:
        output.mkdir(parents=True)
    result={
        "schema":"v130-original-inviscid-SWE-spatial-VTK-recovery-v1",
        "status":"FAIL_OR_UNQUALIFIED_SWE_SPATIAL_RECOVERY",
        "archive":str(archive),"archive_sha256":None,
        "original_native_driver_exitcode":None,
        "original_native_driver_status_preserved":True,
        "original_FE_Hamiltonian_quadratures_recovered":False,
        "original_FE_mass_quadratures_recovered":False,
        "full_time_h_momentum_physical_gate_qualified":False,
        "original_BDF2_nonlinear_solve_modified":False,
        "viscous_rotating_SWE_included":False,
        "errors":[],
    }
    try:
        sha=hashlib.sha256()
        with archive.open("rb") as f:
            for chunk in iter(lambda:f.read(2<<20),b""):
                sha.update(chunk)
        result["archive_sha256"]=sha.hexdigest()
        with ZipFile(archive,"r") as z:
            names=set(z.namelist())
            def read(name):
                if name not in names:
                    raise RuntimeError("Original ZIP lacks "+name)
                return z.read(name)
            manifest=json.loads(read("manifest.json"))
            old=json.loads(read("summary.json"))
            rc=int(read("native.exitcode").strip())
            result["original_native_driver_exitcode"]=rc
            if not (
                manifest["schema"]=="v130-original-SWE-showcase-v1"
                and manifest["source"]==
                    "14181a9431ef621e3df55b7374f61a3746b05e46"
                and manifest["image_revision"]==
                    "1d6d0c90115213298047b28237a7f38087563a35"
                and manifest["model"]=="examples/shallow_water.py"
                and manifest["MPI"]==1
                and manifest["initialization"]=="nonlinear_affine"
                and abs(manifest["final_time"]-.12)<1e-12
                and abs(manifest["dt"]-.0001)<1e-12
                and abs(manifest["save_step"]-.002)<1e-12
                and not manifest["rotating_viscous_research_included"]
                and rc==1
                and old["status"]=="FAIL_OR_UNQUALIFIED_NATIVE_SWE_VISUALS"
            ):
                raise RuntimeError("Archive source/revision/native protocol mismatch")
            log=read("native.host.log").decode("utf-8",errors="replace")
            if not all(m in log for m in (
                "Steps: 1200 (0 rejected, 0 Nonlinear solver failures)",
                "Hamiltonian has been computed",
                "Expected actual linear VTK triangles from GetFEM",
            )):
                raise RuntimeError("Not the native successful-TS/failed-render record")
            with TemporaryDirectory(prefix="scrimp-SWE-recovery-") as tmp:
                base=Path(tmp)
                for field in ("h","p"):
                    prefix=f"native/pv/{field}/"
                    present=sorted(n for n in names if n.startswith(prefix))
                    expected={
                        prefix+field+".pvd",
                        *{prefix+field+f"_{i:05d}.vtu" for i in range(61)},
                    }
                    if set(present)!=expected:
                        raise RuntimeError("Incomplete/extra original GetFEM VTK files "+field)
                    dest=base/"pv"/field
                    dest.mkdir(parents=True)
                    for name in present:
                        path=Path(name)
                        if path.name not in (
                            field+".pvd",
                            *(field+f"_{i:05d}.vtu" for i in range(61)),
                        ):
                            raise RuntimeError("Unsafe original VTK ZIP path")
                        (dest/path.name).write_bytes(read(name))
                hroot=ET.parse(base/"pv/h/h.pvd").getroot()
                original_t=np.array(
                    [float(x.attrib["timestep"])
                     for x in hroot.findall(".//DataSet")],
                    dtype=float,
                )
                if (original_t.shape!=(61,)
                        or not np.all(np.isfinite(original_t))
                        or not np.allclose(
                            original_t,np.arange(61)*.002,
                            atol=1e-7,rtol=0.
                        )):
                    raise RuntimeError("Original saved FE clock incomplete")
                spatial=render_original_shallow_water_spatial_VTU(
                    base,original_t,rho=1000.,
                    max_gif_frames=int(frames),
                    figures_directory=output/"figures",
                )
                if not (
                    spatial["true_FE_states"]==61
                    and spatial["mesh_h_original_quadratic_cells"]==248
                    and spatial["mesh_h_display_linear_subtriangles"]==992
                    and spatial["mesh_p_original_linear_triangles"]==248
                    and len(spatial["original_selected_times"])==
                        spatial["gif_original_FE_frames"]
                ):
                    raise RuntimeError("Incorrect original quadratic GetFEM spatial proof")
                result.update(spatial)
                result["status"]="PASS_RESCUED_NATIVE_SWE_SPATIAL_VTK_ONLY"
                result["source"]=manifest["source"]
                result["original_native_accepted_BDF2_steps"]=1200
                result["original_VTK_h_snapshots"]=61
                result["original_VTK_p_snapshots"]=61
                result["original_saved_final_time"]=float(original_t[-1])
    except Exception as exc:
        result["errors"].append(f"{type(exc).__name__}: {exc}")
    (output/"spatial_recovery.json").write_text(
        json.dumps(result,indent=2)+"\n",encoding="utf-8"
    )
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("archive",type=Path)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--frames",type=int,default=36)
    args=p.parse_args()
    if not 5<=args.frames<=61:
        p.error("frames must be between 5 and 61 real stored FE time instants")
    result=recover(args.archive,args.out,frames=args.frames)
    if result["status"].startswith("PASS_"):
        artifact=Path(str(args.out.resolve())+".zip")
        with ZipFile(artifact,"x",compression=ZIP_DEFLATED,
                     compresslevel=5) as out:
            for f in sorted(args.out.rglob("*")):
                if f.is_file():
                    out.write(f,arcname=str(f.relative_to(args.out)))
        result["recovered_visual_zip"]=str(artifact)
    print(json.dumps(result,indent=2),flush=True)
    raise SystemExit(0 if result["status"].startswith("PASS_") else 1)


if __name__=="__main__":
    main()
