#!/usr/bin/env bash
# SCRIMP v1.3.0: minimal REAL 3D Maxwell with PATCHED Nedelec-II/BDM1.
# Native MUMPS MPI1/2, 10 CN steps, FE energy and original tetra VTK.
# A failed run is always archived and NEVER reported as PASS.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }

docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
  --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
  echo "Incorrect pinned GetFEM/PETSc native image: $revision" >&2; exit 2;
}
python3 -m py_compile examples/maxwell_3d.py scrimp/fem.py tests/test_v130_maxwell3d_source.py tests/test_v130_fem_nedelec_aliases.py
bash -n "$0"
python3 tests/test_v130_fem_nedelec_aliases.py
python3 tests/test_v130_maxwell3d_source.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-maxwell3d-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"v130-maxwell3d-MUMPS-MPI-pilot-v1","source":"%s","image":"%s","image_revision":"%s","canonical":"examples/maxwell_3d.py","domain":"3D tetrahedral Ball, unit radius","fields":"E: Hcurl patched Nedelec-II k1; B: Hdiv patched BDM1","boundary":"natural zero magnetic tangential work (not imposed PEC)","dimensionless_eps_mu":1,"mesh_size":0.7,"final_time":0.2,"time_step":0.02,"steps":10,"not_a_convergence_study":true}\n' \
  "$(git rev-parse HEAD)" "$image" "$revision" > "$out/manifest.json"
echo "MAXWELL 3D: genuine PATCHED FEM_NEDELEC2(3,1) + FEM_BDMK(3,1)"
echo "10 authentic CN steps, MUMPS MPI1/2, H(t)=H(0), tetra VTK E and B"
echo "EVIDENCE: $out"

for ranks in 1 2; do
  name="mumps_mpi$ranks"
  echo "RUN $name"
  mkdir -p "$out/$name" "$out/work_$name"
  rc=0
  "${docker_cmd[@]}" run --rm \
    --user "$(id -u):$(id -g)" \
    --env HOME=/tmp --env MPLBACKEND=Agg \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
    --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
    --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
    --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
    --env SCRIMP_MATRIX_MEMORY_MODE=compact \
    --env SCRIMP_MUMPS_RESET_POLICY=on_change \
    --env SCRIMP_GETFEM_RHS_MODE=auto \
    --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
    --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
    --volume "$root/examples/maxwell_3d.py:/workspace/scrimp/examples/maxwell_3d.py:ro" \
    --volume "$out/$name:/artifacts:rw" \
    --volume "$out/work_$name:/scrimp-work:rw" \
    --workdir /scrimp-work \
    "$image" timeout -k 5s "${SCRIMP_V130_MAXWELL3D_TIMEOUT:-360s}" \
    mpiexec --oversubscribe --bind-to none -n "$ranks" \
    python -u -m examples.maxwell_3d \
      --solver mumps --mesh-size 0.7 --final-time 0.2 --time-step 0.02 \
      --check --export --output /artifacts \
    >"$out/$name.host.log" 2>&1 || rc=$?
  echo "$rc" >"$out/$name.exitcode"
  if [[ "$rc" != 0 ]]; then
    echo "UNQUALIFIED $name exit=$rc; retaining native error evidence" >&2
    tail -n 45 "$out/$name.host.log" >&2
  fi
done

python3 - "$out" <<'PY'
import base64
import json
import struct
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
import numpy as np

folder = Path(sys.argv[1])
names = ("mumps_mpi1", "mumps_mpi2")
keys = ("t", "z", "H", "electric_energy", "magnetic_energy",
        "div_B_squared", "step_energy_defect",
        "cumulative_energy_defect")
summary = {
    "schema": "v130-maxwell3d-Nedelec-RT0-native-v1",
    "status": "FAIL_OR_UNQUALIFIED_MAXWELL3D_NATIVE",
    "cases": {}, "MPI1_vs_MPI2": {}, "errors": [],
    "superlu_dist_qualified": False,
    "mesh_time_convergence_qualified": False,
    "published_sources_unchanged": True,
}
for ranks,name in enumerate(names, start=1):
    summary["cases"][name] = {
        "MPI": ranks,
        "exit_code": int((folder/f"{name}.exitcode").read_text())
    }
try:
    case_traces = {}
    for ranks,name in enumerate(names, start=1):
        rec=summary["cases"][name]
        if rec["exit_code"]:
            raise RuntimeError(f"{name}: native exit {rec['exit_code']}")
        rank_traces = []
        for rank in range(ranks):
            case_dir = folder/name
            meta=json.loads((case_dir/f"rank{rank}.json").read_text())
            with np.load(case_dir/f"rank{rank}.npz",allow_pickle=False) as data:
                sample={k:data[k].copy() for k in data.files}
            gates=meta.get("physical_gates",{})
            # Ensure this is the PATCHED 3D FEM pair, not the previously
            # accepted baseline GetFEM Nedelec0/RT0 result.
            families=meta.get("FE",{})
            expected_fe={
                "E": ("NED2", 1, "FEM_NEDELEC2(3,1)"),
                "B": ("BDM", 1, "FEM_BDMK(3,1)"),
            }
            if set(families) != set(expected_fe):
                raise RuntimeError(f"{name}/rank{rank}: missing patched FEM metadata")
            for label, expected in expected_fe.items():
                item=families[label]
                if ((item.get("family"), item.get("k"),item.get("getfem"))
                    != expected or int(item.get("DOFs",0)) <= 0):
                    raise RuntimeError(f"{name}/rank{rank}: wrong patched {label} FE")
            if sum(families[f]["DOFs"] for f in expected_fe) != meta.get("DOFs"):
                raise RuntimeError(f"{name}/rank{rank}: FE DOFs do not add up")
            if (meta.get("status") != "PASS_NATIVE_MAXWELL3D"
                or meta.get("MPI") != ranks or meta.get("rank") != rank
                or meta.get("solver") != "mumps"
                or meta.get("steps") != 10 or meta.get("TS_reason",0) <= 0
                or len(gates) != 9 or not all(gates.values())
                or any(k not in sample for k in keys)):
                raise RuntimeError(f"{name}/rank{rank}: native physics gate FAIL")
            t=sample["t"];H=sample["H"];z=sample["z"]
            if (t.shape!=(11,) or H.shape!=(11,)
                or z.shape!=(11,meta["DOFs"])
                or not np.allclose(t,np.arange(11)*.02,rtol=0,atol=2e-10)
                or any(not np.all(np.isfinite(x)) for x in sample.values())):
                raise RuntimeError(f"{name}/rank{rank}: invalid FE history")
            if (np.max(np.abs(H-H[0])) > 1e-9
                or np.max(np.abs(H-sample["electric_energy"]-
                                 sample["magnetic_energy"])) > 1e-9
                or sample["magnetic_energy"][-1] <= 1e-6*H[0]
                or np.max(sample["div_B_squared"]) > 1e-8*max(1.,float(H[0]))
                or not np.allclose(sample["step_energy_defect"],
                                   np.diff(H),atol=1e-10,rtol=0)):
                raise RuntimeError(f"{name}/rank{rank}: independent energy FAIL")
            rank_traces.append(sample)
        for other in rank_traces[1:]:
            for key in keys:
                if not np.allclose(other[key],rank_traces[0][key],
                                   atol=5e-9,rtol=2e-8):
                    raise RuntimeError(f"{name}: rank replica mismatch {key}")

        for field in ("E","B"):
            pvd=folder/name/"pv"/field/f"{field}.pvd"
            frames=ET.parse(pvd).findall(".//DataSet")
            if len(frames)!=11:
                raise RuntimeError(f"{name}: {field} PVD needs 11 native frames")
            for i,item in enumerate(frames):
                if abs(float(item.get("timestep"))-rank_traces[0]["t"][i])>1e-10:
                    raise RuntimeError(f"{name}: {field} invalid native time {i}")
                vtuname=Path(item.attrib["file"])
                if vtuname.name!=str(vtuname):
                    raise RuntimeError(f"{name}: unsafe VTU relative name")
                root=ET.fromstring((pvd.parent/vtuname).read_bytes())
                piece=root.find(".//Piece")
                if (piece is None or int(piece.attrib["NumberOfPoints"]) <= 0
                    or int(piece.attrib["NumberOfCells"]) <= 0
                    or piece.find(f"./PointData/DataArray[@Name='{field}']")
                       is None):
                    raise RuntimeError(f"{name}: invalid 3D {field} VTU/{i}")
                # REAL VTK tetrahedra, not just a file whose mesh.dim()=3:
                # GetFEM binary cell types use an int32 byte-size prefix.
                cell_node=piece.find("./Cells/DataArray[@Name='types']")
                if (cell_node is None or cell_node.get("type") != "Int32"
                    or cell_node.get("format") != "binary"):
                    raise RuntimeError(f"{name}: missing binary 3D cell types")
                payload=base64.b64decode((cell_node.text or "").strip())
                if len(payload)<4 or struct.unpack("<I",payload[:4])[0]!=len(payload)-4:
                    raise RuntimeError(f"{name}: corrupt native 3D VTK cells")
                cell_types=np.frombuffer(payload[4:],dtype="<i4")
                if (len(cell_types)!=int(piece.attrib["NumberOfCells"])
                    or not np.all(np.isin(cell_types, (10,24)))):
                    raise RuntimeError(
                        f"{name}: {field} exported non-tetrahedral GetFEM cells"
                    )
        # --export now calls the renderer INSIDE examples/maxwell_3d.py.
        # Strictly accept only real files and the provenance of this very
        # physical FE trajectory; no post-hoc external renderer is used.
        visuals=folder/name/"visuals"
        prov=json.loads((visuals/"provenance.json").read_text())
        expected_assets=(
            "maxwell3d_echange_energie.png",
            "maxwell3d_echange_energie.svg",
            "maxwell3d_vecteurs_E_tfinal.png",
            "maxwell3d_vecteurs_B_tfinal.png",
            "maxwell3d_B_3D_animation.gif",
        )
        rank0_meta=json.loads((folder/name/"rank0.json").read_text())
        if (prov.get("schema")!="scrimp-maxwell3d-native-visuals-v1"
            or prov.get("MPI_writer_rank")!=0
            or prov.get("frames")!=11 or prov.get("tetrahedra")!=110
            or prov.get("points")!=58
            or prov.get("original_FE_state_SHA256")!=rank0_meta.get("state_SHA256")
            or tuple(prov.get("assets",()))!=expected_assets
            or rank0_meta.get("visuals")!=prov):
            raise RuntimeError(f"{name}: visual provenance does not match real FE")
        for asset in expected_assets:
            path=visuals/asset
            if not path.is_file() or path.stat().st_size<128:
                raise RuntimeError(f"{name}: missing/empty native visual {asset}")
            signature=path.read_bytes()[:12]
            if asset.endswith(".png") and not signature.startswith(b"\x89PNG\r\n\x1a\n"):
                raise RuntimeError(f"{name}: invalid PNG file {asset}")
            if asset.endswith(".gif") and signature[:6] not in (b"GIF87a",b"GIF89a"):
                raise RuntimeError(f"{name}: invalid GIF file {asset}")
            if asset.endswith(".svg"):
                ET.parse(path)
        # Pillow is present with Matplotlib in the pinned Docker image;
        # use it to verify the GIF includes *all* 11 native field states.
        from PIL import Image
        with Image.open(visuals/expected_assets[-1]) as gif:
            if gif.n_frames!=11:
                raise RuntimeError(f"{name}: animated GIF lost native CN frames")
        if not (visuals/"README.txt").is_file():
            raise RuntimeError(f"{name}: missing in-example figure provenance notes")
        rec.update({"status":"PASS_NATIVE_MAXWELL3D_MUMPS_CN_VTK_VISUALS",
                    "DOFs":int(rank_traces[0]["z"].shape[1]),
                    "steps":10,"VTK_fields":{"E":11,"B":11},
                    "visuals":list(expected_assets),
                    "visuals_generated_by":"examples/maxwell_3d.py",
                    "max_energy_defect":float(np.max(np.abs(
                        rank_traces[0]["cumulative_energy_defect"]))),
                    "rank_traces_agree":True})
        case_traces[name] = rank_traces[0]
    a,b=(case_traces[name] for name in names)
    for k in keys:
        if a[k].shape!=b[k].shape:
            raise RuntimeError(f"MPI1/2 different FE dimensions for {k}")
        d=float(np.max(np.abs(a[k]-b[k])))
        summary["MPI1_vs_MPI2"][k]=d
        if not np.allclose(a[k],b[k],atol=5e-9,rtol=2e-8):
            raise RuntimeError(f"MPI1/2 FE difference too large for {k}: {d}")
    summary["status"]="PASS_MAXWELL3D_MUMPS_MPI1_MPI2"
except Exception as exc:
    summary["errors"].append(f"{type(exc).__name__}: {exc}")

(folder/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
archive=Path(str(folder)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as out:
    for file in sorted(folder.rglob("*")):
        if file.is_file():
            out.write(file,arcname=str(file.relative_to(folder)))
print(json.dumps({**summary,"archive":str(archive),
                  "evidence":str(folder)},indent=2),flush=True)
raise SystemExit(0 if summary["status"].startswith("PASS_") else 1)
PY
