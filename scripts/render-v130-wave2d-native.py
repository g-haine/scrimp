#!/usr/bin/env python3
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
"""User-facing images/animation from ORIGINAL GetFEM VTK and CN power archives.

Supports 2A, 2B, 2C when they provide the native rank0.npz / pv/*.pvd
evidence layout. Field images use original VTK Float32 samples for
visualization ONLY; passivity and energy audits use raw GetFEM NPZ.
No reconstructed displacement, smoothed time history or synthetic fields.
"""
import argparse
import base64
import csv
import hashlib
import json
import struct
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np


class Evidence:
    def __init__(self, path):
        self.path = Path(path)
        self.zip = zipfile.ZipFile(self.path) if self.path.is_file() else None
        if self.zip is None and not self.path.is_dir():
            raise FileNotFoundError(self.path)

    def read(self, relative):
        relative = str(relative).replace("\\", "/")
        if self.zip is not None:
            return self.zip.read(relative)
        return (self.path / relative).read_bytes()

    def exists(self, relative):
        try:
            self.read(relative)
            return True
        except (FileNotFoundError, KeyError):
            return False

    def json(self, relative):
        return json.loads(self.read(relative).decode("utf-8"))


def vtk_binary(node):
    if node is None or node.get("format") != "binary":
        raise ValueError("Expected original GetFEM inline-binary VTK DataArray")
    dtype = {
        "Float32": "<f4", "Float64": "<f8",
        "Int32": "<i4", "Int64": "<i8", "UInt8": "u1",
    }.get(node.get("type"))
    if dtype is None:
        raise ValueError("Unsupported VTK numeric dtype")
    raw = base64.b64decode((node.text or "").strip())
    if len(raw) < 4:
        raise ValueError("Truncated binary VTK header")
    nbytes = struct.unpack("<I", raw[:4])[0]
    if len(raw) != nbytes + 4:
        raise ValueError("GetFEM VTK byte-count mismatch")
    out = np.frombuffer(raw, dtype=np.dtype(dtype), offset=4)
    components = int(node.get("NumberOfComponents", "1"))
    if components > 1:
        if len(out) % components:
            raise ValueError("Invalid VTK component count")
        out = out.reshape((-1, components))
    return out.copy()


def vtk_snapshot(evidence, case, field, filename):
    file = f"{case}/pv/{field}/{filename}"
    node = ET.fromstring(evidence.read(file))
    points = vtk_binary(node.find(".//Points/DataArray"))
    connectivity = vtk_binary(node.find(".//Cells/DataArray[@Name='connectivity']"))
    offsets = vtk_binary(node.find(".//Cells/DataArray[@Name='offsets']"))
    types = vtk_binary(node.find(".//Cells/DataArray[@Name='types']"))
    values = vtk_binary(node.find(f".//PointData/DataArray[@Name='{field}']"))
    if points.ndim != 2 or points.shape[1] < 2:
        raise ValueError("Unexpected GetFEM VTK coordinates")
    if len(types) != len(offsets):
        raise ValueError("Invalid VTK cell offsets")
    if np.any(~np.isin(types, (5, 22))):
        raise ValueError("Only native GetFEM triangular cells supported")
    starts = np.r_[0, offsets[:-1]]
    triangle = np.array([connectivity[a:b][:3]
                         for a, b in zip(starts, offsets)], dtype=int)
    if triangle.ndim != 2 or triangle.shape[1] != 3:
        raise ValueError("Invalid triangle connectivity")
    if np.any(triangle < 0) or np.any(triangle >= len(points)):
        raise ValueError("Invalid VTK vertex index")
    if len(values) != len(points):
        raise ValueError("VTK FE field / mesh size mismatch")
    if values.ndim == 2:
        scalar = np.linalg.norm(values, axis=1)
    else:
        scalar = values
    if not np.all(np.isfinite(scalar)):
        raise ValueError("Nonfinite authentic VTK field")
    return points[:, :2], triangle, scalar


def field_history(evidence, case, field):
    root = ET.fromstring(evidence.read(f"{case}/pv/{field}/{field}.pvd"))
    frames = root.findall(".//DataSet")
    if len(frames) < 2:
        raise ValueError(f"No time series for {field}")
    result = []
    for frame in frames:
        name = Path(frame.attrib["file"])
        if name.name != str(name):
            raise ValueError("Unsafe PVD relative file path")
        t = float(frame.attrib["timestep"])
        xy, triangles, values = vtk_snapshot(evidence, case, field, name.name)
        result.append((t, xy, triangles, values))
    if not np.all(np.diff([r[0] for r in result]) > 0):
        raise ValueError("Unordered native VTK timestamps")
    return result


def normalized_balance(raw):
    t = raw["t"].astype(float)
    H = raw["H"].astype(float)
    if "cumulative_internal_loss" in raw:
        # 2C: closed system, external work = 0; boundary feedback dissipates.
        W = np.zeros_like(H)
        Li = raw["cumulative_internal_loss"]
        Lb = raw["cumulative_feedback_loss"]
        pb = raw["power_boundary"]
        pi = raw["dissipation_internal"]
        kind = "2C · Coenergy: inner disk and boundary feedback"
    else:
        W = raw["signed_work"] if "signed_work" in raw else np.zeros_like(H)
        Li = raw["cumulative_loss"] if "cumulative_loss" in raw else np.zeros_like(H)
        Lb = np.zeros_like(H)
        pb = raw["power_left"] if "power_left" in raw else np.zeros_like(H)
        pi = raw["dissipated_power"] if "dissipated_power" in raw else np.zeros_like(H)
        kind = ("2B · Dissipative membrane" if "cumulative_loss" in raw
                else "2A · Conservative membrane")
    if any(np.shape(a) != np.shape(t) for a in (H, W, Li, Lb, pb, pi)):
        raise ValueError("Malformed genuine CN energy/power observations")
    if not all(np.all(np.isfinite(a)) for a in (t,H,W,Li,Lb,pb,pi)):
        raise ValueError("Nonfinite original energy data")
    if len(t) < 2 or not np.all(np.diff(t)>0):
        raise ValueError("Invalid accepted CN times")
    return dict(t=t,H=H,W=W,Li=Li,Lb=Lb,Pboundary=pb,Pinternal=pi,
                residual=H-H[0]-W+Li+Lb,kind=kind)


def energy_figures(out, trace):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t=trace["t"];H=trace["H"];W=trace["W"];Li=trace["Li"];Lb=trace["Lb"]
    fig, axes=plt.subplots(3,1,figsize=(9,9),sharex=True,layout="constrained")
    axes[0].plot(t,H-H[0],"o-",label="H(t) - H(0) · genuine GetFEM")
    axes[0].plot(t,W-Li-Lb,"--",label="W external − losses")
    axes[0].plot(t,Li,":",label="Internal dissipated energy")
    if np.any(Lb):
        axes[0].plot(t,Lb,":",label="Boundary feedback dissipated energy")
    axes[0].set_ylabel("Energy")
    axes[0].legend()
    axes[1].plot(t,trace["Pboundary"],"o-",label="Signed boundary power")
    axes[1].plot(t,trace["Pinternal"],"-",label="Positive inner/viscous dissipation")
    axes[1].set_ylabel("Power")
    axes[1].legend()
    axes[2].plot(t,trace["residual"],"o-")
    axes[2].set_ylabel("H − H0 − W + losses")
    axes[2].set_xlabel("Time")
    for ax in axes:ax.grid(True,alpha=.25)
    fig.suptitle(trace["kind"]+" · authentic CN power/energy observations")
    for ext in (".png",".svg"):
        fig.savefig(out/("energy_balance"+ext),dpi=160)
    plt.close(fig)


def render_maps(out,evidence,case,fields,trace,no_gif=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib.animation import FuncAnimation, PillowWriter
    histories={f:field_history(evidence,case,f) for f in fields}
    for field,frames in histories.items():
        if len(frames)!=len(trace["t"]) or any(
            abs(frame[0]-trace["t"][i])>1e-10
            for i,frame in enumerate(frames)):
            raise ValueError(f"Native {field} VTK and CN history out of sync")
    times=(0,len(trace["t"])//2,len(trace["t"])-1)
    fig, axes=plt.subplots(len(fields),3,figsize=(12,3.7*len(fields)),
                           squeeze=False,layout="constrained")
    for row, field in enumerate(fields):
        frames=histories[field]
        scalar_faces=[np.mean(x[3][x[2]],axis=1) for x in frames]
        lo=min(float(np.min(s)) for s in scalar_faces)
        hi=max(float(np.max(s)) for s in scalar_faces)
        if hi<=lo:hi=lo+1e-12
        for col,index in enumerate(times):
            t,xy,tri,_=frames[index]
            artists=PolyCollection(xy[tri],array=scalar_faces[index],
                                   cmap="viridis",edgecolors="none")
            artists.set_clim(lo,hi)
            ax=axes[row,col]
            ax.add_collection(artists)
            ax.autoscale()
            ax.set_aspect("equal",adjustable="box")
            ax.set_title(f"{field}"+(" · magnitude" if frames[index][3].ndim==2 else "")+
                         f"  · t={t:.3f}")
            ax.set_xlabel("x");ax.set_ylabel("y")
            fig.colorbar(artists,ax=ax,shrink=.75)
    fig.suptitle(trace["kind"]+" · original GetFEM VTK fields (Float32)")
    for ext in (".png",".svg"):
        fig.savefig(out/("field_snapshots"+ext),dpi=160)
    plt.close(fig)

    if not no_gif and "p" in histories:
        frames=histories["p"]
        faces=[np.mean(x[3][x[2]],axis=1) for x in frames]
        lo=min(float(np.min(s)) for s in faces)
        hi=max(float(np.max(s)) for s in faces)
        if hi<=lo:hi=lo+1e-12
        fig,ax=plt.subplots(figsize=(6,6),layout="constrained")
        mesh=PolyCollection(frames[0][1][frames[0][2]],array=faces[0],
                            cmap="viridis",edgecolors="none")
        mesh.set_clim(lo,hi)
        ax.add_collection(mesh);ax.autoscale();ax.set_aspect("equal")
        ax.set_xlabel("x");ax.set_ylabel("y")
        fig.colorbar(mesh,ax=ax,label="p · original FE VTK sample")
        title=ax.set_title("")
        def draw(k):
            mesh.set_verts(frames[k][1][frames[k][2]])
            mesh.set_array(faces[k])
            title.set_text(f"{trace['kind']}  |  t={frames[k][0]:.3f}")
            return mesh,title
        movie=FuncAnimation(fig,draw,frames=len(frames),interval=350,blit=False)
        movie.save(out/"p_evolution.gif",writer=PillowWriter(fps=1000/350),
                   dpi=110)
        plt.close(fig)
    return {k:len(v) for k,v in histories.items()}


def write_csv(out,trace):
    with (out/"energy_and_ports.csv").open("w",newline="") as file:
        writer=csv.writer(file)
        writer.writerow(("t","H","H_minus_H0","external_signed_work",
                         "internal_cumulative_loss","feedback_cumulative_loss",
                         "signed_boundary_power","internal_dissipative_power",
                         "energy_balance_defect"))
        for i,t in enumerate(trace["t"]):
            writer.writerow((float(t),float(trace["H"][i]),
                            float(trace["H"][i]-trace["H"][0]),
                            float(trace["W"][i]),float(trace["Li"][i]),
                            float(trace["Lb"][i]),float(trace["Pboundary"][i]),
                            float(trace["Pinternal"][i]),
                            float(trace["residual"][i])))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,required=True,
                        help="Native pilot archive ZIP or original output directory")
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--case",default="mumps_mpi1")
    parser.add_argument("--no-gif",action="store_true")
    args=parser.parse_args()
    data=Evidence(args.input)
    if not data.exists("manifest.json"):
        raise ValueError("Missing authentic native experiment provenance")
    manifest=data.json("manifest.json")
    if data.exists("summary.json"):
        summary=data.json("summary.json")
        if not str(summary.get("status","")).startswith("PASS_"):
            raise ValueError("Refusing to render unqualified native physics")
        if (summary.get("cases",{}).get(args.case,{}).get("exit_code")!=0):
            raise ValueError("Selected native case failed")
    if data.exists(f"{args.case}/rank0.json"):
        meta=data.json(f"{args.case}/rank0.json")
        if not str(meta.get("status","")).startswith("PASS_"):
            raise ValueError("Selected rank0 original native physical gates failed")
    import io
    with np.load(io.BytesIO(data.read(f"{args.case}/rank0.npz")),
                 allow_pickle=False) as loaded:
        raw={name:loaded[name] for name in loaded.files}
    trace=normalized_balance(raw)
    fields=[f for f in ("p","q","e_r") if data.exists(
        f"{args.case}/pv/{f}/{f}.pvd")]
    if "p" not in fields or "q" not in fields:
        raise ValueError("Authentic q and p VTK fields required for visual outputs")
    args.output.mkdir(parents=True,exist_ok=True)
    write_csv(args.output,trace)
    energy_figures(args.output,trace)
    vtk=render_maps(args.output,data,args.case,fields,trace,args.no_gif)
    provenance={
        "status":"RENDERED_FROM_REAL_NATIVE_GETFEM_FE_AND_VTK",
        "source":manifest.get("source"),
        "image_revision":manifest.get("image_revision"),
        "case":args.case,"model":trace["kind"],
        "VTK_fields":vtk,"frames":len(trace["t"]),
        "input_SHA256":(hashlib.sha256(args.input.read_bytes()).hexdigest()
                        if args.input.is_file() else None),
        "visual_field_precision":"GetFEM VTU Float32; visual only",
        "energy_precision":"Original FE GetFEM Hamiltonian and CN port powers",
        "not_an_independent_physical_qualification":True,
    }
    (args.output/"provenance.json").write_text(
        json.dumps(provenance,indent=2)+"\n",encoding="utf-8")
    (args.output/"README.txt").write_text(
        "SCRIMP wave 2D: genuine accepted VTK/FE outputs.\n"
        "Open the native pv/<field>/<field>.pvd time series in ParaView.\n"
        "field_snapshots.png/.svg: native p, q (vector magnitude), e_r when present.\n"
        "p_evolution.gif: true p snapshots at native saved times.\n"
        "energy_balance.png/.svg and energy_and_ports.csv: original CN powers.\n"
        "Field VTUs are Float32, not scientific energy quadrature.\n"
        "A short pilot demonstrates the pipeline, not visually developed waves.\n",
        encoding="utf-8")
    print(json.dumps({"status":provenance["status"],
                      "output":str(args.output),
                      "assets":sorted(p.name for p in args.output.iterdir())},
                     indent=2),flush=True)


if __name__=="__main__":
    main()
