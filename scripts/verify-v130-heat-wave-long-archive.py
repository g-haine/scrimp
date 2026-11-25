#!/usr/bin/env python3
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
"""Re-audit ORIGINAL long Heat-Wave native ZIP evidence, without rerunning CN.

The 2026-11-25 b9d56fcc pilot finished both genuine 15-second simulations
(exit code 0) but its POST-RUN PVD check used XML_Element["file"] rather
than XML_Element.attrib["file"]. Verify the immutable existing artifacts
rather than waste CPU recomputing a completed 3000-step GetFEM run.

The result is a NEW sidecar JSON report, NOT a rewrite of the old ZIP or
of its historical FAIL status. A visualization PASS does NOT certify the
CN heat+wave interface power balance at the full 15-second horizon,
nor a whispering-gallery eigenmode.
"""
import argparse
import hashlib
from io import BytesIO
import json
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

import numpy as np
from PIL import Image


def check(archive):
    archive=Path(archive)
    out={
        "schema":"v130-heat-wave-long-archive-reaudit-v1",
        "status":"FAIL_OR_UNQUALIFIED_ARCHIVED_LONG_VISUALS",
        "original_archive":str(archive),
        "archive_sha256":None,
        "original_native_exit_code":None,
        "profile":None,
        "full_horizon_energy_power_identity_qualified":False,
        "whispering_gallery_eigenmode_qualified":False,
        "errors":[],
    }
    try:
        if not archive.is_file():
            raise FileNotFoundError(archive)
        sha=hashlib.sha256()
        with archive.open("rb") as fh:
            while chunk:=fh.read(8*1024*1024):
                sha.update(chunk)
        out["archive_sha256"]=sha.hexdigest()
        with ZipFile(archive,"r") as zip:
            names=set(zip.namelist())
            def item(name, *, minimum=1):
                if name not in names or zip.getinfo(name).file_size<minimum:
                    raise RuntimeError("Missing/empty authentic ZIP member: "+name)
                return zip.read(name)
            manifest=json.loads(item("manifest.json"))
            previous=json.loads(item("summary.json"))
            exit_code=int(item("native.exitcode").decode().strip())
            out["original_native_exit_code"]=exit_code
            profile=manifest["profile"]
            out["profile"]=profile
            if not (
                manifest["schema"]=="v130-heat-wave-long-showcase-v1"
                and manifest["source"]=="b9d56fccc18c1f9facacdd9b6c2645c4de25fc47"
                and manifest["revision"]=="1d6d0c90115213298047b28237a7f38087563a35"
                and manifest["MPI"]==1 and manifest["time_scheme"]=="CN"
                and manifest["dt"]==0.005 and manifest["tf"]==15.
                and manifest["dt_save"]==0.05
                and profile in ("gaussian","gallery_packet")
                and previous["native_exit_code"]==exit_code==0
                and previous["profile"]==profile
                and previous["errors"]==["TypeError: element indices must be integers"]
            ):
                raise RuntimeError("Original native source/image or failure diagnosis differs")
            prefix="work/outputs/"
            root=prefix+"showcase/heat_wave_annulus/"
            evidence=json.loads(item(root+"provenance.json"))
            with np.load(BytesIO(item(root+"heat_wave_FE_H_long.npz")),
                         allow_pickle=False) as npz:
                t=np.asarray(npz["t"],dtype=np.float64)
                H=np.asarray(npz["H"],dtype=np.float64)
            if not (
                evidence["schema"]=="scrimp-heat-wave-long-native-showcase-v1"
                and evidence["initial_profile"]==profile
                and evidence["no_time_interpolation"] is True
                and evidence["whispering_gallery_eigenmode_certified"] is False
                and evidence["full_field_energy_localization_certified"] is False
                and t.shape==H.shape==(301,)
                and np.all(np.isfinite(H)) and np.min(H)>0.
                and np.allclose(t,np.arange(301)*0.05,rtol=0,atol=1e-9)
                and abs(H[0]-evidence["initial_H"])<1e-12
                and abs(H[-1]-evidence["final_H"])<1e-12
            ):
                raise RuntimeError("GetFEM time/H data not original or inconsistent")
            max_increment=float(np.max(np.diff(H)))
            if max_increment>1e-7*max(1.,float(H[0])):
                raise RuntimeError("Long GetFEM Hamiltonian increases beyond audit tolerance")
            graphics=[
                "heat_wave_H_long_loglog.png",
                "heat_wave_H_long_loglog.svg",
                "wave_annulus_long.gif",
                "wave_outer_shell_kinetic_proxy.png",
            ]
            for filename in graphics:
                raw=item(root+filename,minimum=250)
                if filename.endswith(".png") and raw[:8]!=b"\x89PNG\r\n\x1a\n":
                    raise RuntimeError("Invalid native PNG: "+filename)
                if filename.endswith(".gif") and raw[:6] not in (b"GIF89a",b"GIF87a"):
                    raise RuntimeError("Invalid native GIF: "+filename)
                if filename.endswith(".svg"):
                    ET.fromstring(raw)
            # Four genuine selected VTK snapshots are also required.
            snapshots=[
                n for n in names
                if n.startswith(root+"wave_annulus_t") and n.endswith(".png")
            ]
            if len(snapshots)!=4:
                raise RuntimeError("Expected four native VTK wave snapshots")
            for filename in snapshots:
                if item(filename,minimum=250)[:8]!=b"\x89PNG\r\n\x1a\n":
                    raise RuntimeError("Invalid native VTK wave PNG "+filename)
            gif=Image.open(BytesIO(item(root+"wave_annulus_long.gif")))
            nframes=gif.n_frames
            if nframes!=48 or evidence["gif_frames_real"]!=nframes:
                raise RuntimeError("GIF frame count differs from saved FE evidence")
            selected=np.asarray(evidence["selected_true_FE_times"],dtype=float)
            expected_indices=np.unique(np.linspace(0,300,48,dtype=int))
            if (selected.shape!=(48,)
                    or not np.allclose(selected,t[expected_indices],rtol=0,atol=1e-9)):
                raise RuntimeError("GIF timestamps not true selected GetFEM states")
            concentration=evidence["kinetic_proxy_shell_fraction"]
            if len(concentration)!=nframes:
                raise RuntimeError("Kinetic proxy trace does not match GIF clock")
            fractions=np.array(
                [np.nan if x is None else float(x) for x in concentration],
                dtype=float,
            )
            if (np.any(np.isinf(fractions))
                    or np.any(fractions[np.isfinite(fractions)]<0.)
                    or np.any(fractions[np.isfinite(fractions)]>1.)):
                raise RuntimeError("Invalid bounded annular kinetic proxy fraction")

            counts={}
            for field in ("p","q","T"):
                pvd_dir=prefix+f"pv/{field}/"
                pvd=ET.fromstring(item(pvd_dir+field+".pvd"))
                frames=pvd.findall(".//DataSet")
                if len(frames)!=len(t):
                    raise RuntimeError("GetFEM PVD count differs from true FE clock: "+field)
                for i,frame in enumerate(frames):
                    # Correct API: Element.attrib, NOT Element["file"].
                    rel=Path(frame.attrib["file"])
                    if rel.name!=str(rel) or rel.suffix!=".vtu":
                        raise RuntimeError("Illegal VTK PVD member path")
                    if abs(float(frame.attrib["timestep"])-float(t[i]))>1e-9:
                        raise RuntimeError("VTK PVD/FE clock mismatch: "+field)
                    # Use only central ZIP directory metadata for all
                    # 903 native frame file existence checks.
                    name=pvd_dir+str(rel)
                    if name not in names or zip.getinfo(name).file_size<128:
                        raise RuntimeError("Missing actual GetFEM VTU: "+name)
                for i in (0,len(frames)//2,len(frames)-1):
                    rel=frames[i].attrib["file"]
                    tree=ET.fromstring(item(pvd_dir+rel))
                    piece=tree.find(".//Piece")
                    node=piece.find(f"./PointData/DataArray[@Name='{field}']") if piece is not None else None
                    if (piece is None or node is None
                            or int(piece.attrib["NumberOfCells"])<=0
                            or int(piece.attrib["NumberOfPoints"])<=0
                            or node.get("format") not in ("ascii","binary")
                            or not (node.text or "").strip()):
                        raise RuntimeError("Invalid genuine GetFEM VTU field: "+field)
                counts[field]=len(frames)
            out.update({
                "status":"PASS_ARCHIVED_LONG_NATIVE_VISUALS_NOT_WG_MODE_PROOF",
                "times_saved":len(t),
                "vtus_per_field":counts,
                "gif_frames_real":nframes,
                "H_initial":float(H[0]),
                "H_final":float(H[-1]),
                "relative_H_loss":float(1-H[-1]/H[0]),
                "max_H_increment":max_increment,
                "outer_shell_kinetic_fraction_mean":float(np.nanmean(fractions)),
                "outer_shell_kinetic_fraction_final":float(fractions[-1]),
                "original_archive_status_preserved":previous["status"],
                "graphics_from_canonical_python":True,
            })
    except Exception as exc:
        out["errors"].append(f"{type(exc).__name__}: {exc}")
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archives",type=Path,nargs="+")
    args=parser.parse_args()
    passed=True
    for archive in args.archives:
        result=check(archive)
        # Sidecar for the *new* independent audit; never alter the ZIP,
        # its historical FAIL or the original scientific evidence.
        sidecar=Path(str(archive)+".reaudit.json")
        sidecar.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
        print(json.dumps({**result,"reaudit_file":str(sidecar)},indent=2))
        passed &= result["status"].startswith("PASS_")
    raise SystemExit(0 if passed else 1)


if __name__=="__main__":
    main()
