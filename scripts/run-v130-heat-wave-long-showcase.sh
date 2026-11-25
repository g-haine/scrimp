#!/usr/bin/env bash
# v1.3.0: REAL long-time Heat–Wave example and annular-wave illustrations.
# This is a separate observational SHOWCASE, not the strict CN physics gate.
# It never reports a whispering-gallery eigenmode without modal analysis.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }
docker_cmd=(sudo docker)
[[ "${SCRIMP_DOCKER_NO_SUDO:-0}" != 1 ]] || docker_cmd=(docker)
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
  --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == 1d6d0c90115213298047b28237a7f38087563a35 ]] || {
  echo "Wrong pinned scientific image $revision" >&2; exit 2;
}
profile="${SCRIMP_V130_HW_SHOWCASE_PROFILE:-gallery_packet}"
case "$profile" in
  gallery_packet|gaussian) ;;
  *) echo "Profile must be gallery_packet or gaussian" >&2; exit 2;;
esac
python3 -m py_compile examples/heat_wave.py scrimp/utils/fe_initialization.py \
  tests/test_v130_heat_wave_showcase_source.py
python3 tests/test_v130_heat_wave_showcase_source.py
python3 scripts/check-v130-published-immutable.py
bash -n "$0"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-heat-wave-long-$profile-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
mkdir "$out/work"
printf '{"schema":"v130-heat-wave-long-showcase-v1","source":"%s","image":"%s","revision":"%s","profile":"%s","waves":"outer_annulus","heat":"central_disk","time_scheme":"CN","dt":0.005,"tf":15.0,"dt_save":0.05,"MPI":1,"is_modal_eigenmode_proof":false,"physical_short_gate_separate":true}\n' \
  "$(git rev-parse HEAD)" "$image" "$revision" "$profile" > "$out/manifest.json"
echo "Heat–Wave LONG native: t_f=15 dt=.005 saved .05; profile=$profile"
echo "Native log–log H + actual annular p/q GetFEM PVD + p GIF + shell diagnostic"
echo "WHISPERING GALLERY IS A CANDIDATE VISUAL PATTERN, NOT A CERTIFIED MODE"
echo "EVIDENCE: $out"
rc=0
"${docker_cmd[@]}" run --rm \
  --user "$(id -u):$(id -g)" --env HOME=/tmp \
  --env MPLBACKEND=Agg --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 \
  --env PYTHONPATH=/workspace/scrimp \
  --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
  --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
  --env SCRIMP_MATRIX_MEMORY_MODE=compact \
  --env SCRIMP_MUMPS_RESET_POLICY=on_change \
  --env SCRIMP_GETFEM_RHS_MODE=auto \
  --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
  --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
  --volume "$root/examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro" \
  --volume "$out/work:/scrimp-work:rw" \
  --workdir /scrimp-work \
  "$image" timeout -k 10s "${SCRIMP_V130_HW_LONG_TIMEOUT:-7200s}" \
  mpiexec --oversubscribe --bind-to none -n 1 \
  python -u -m examples.heat_wave \
    --heat-region 1 --wave-region 2 \
    --initial-wave-boundary clamp \
    --initial-wave-profile "$profile" \
    --final-time 15. --time-step .005 --save-step .05 \
    --factor-solver mumps --long-showcase \
  > "$out/native.host.log" 2>&1 || rc=$?
echo "$rc" > "$out/native.exitcode"
if [[ "$rc" != 0 ]]; then
  echo "UNQUALIFIED long native exit=$rc; archiving partial evidence" >&2
  tail -n 65 "$out/native.host.log" >&2
fi

python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
import numpy as np
from PIL import Image

root=Path(sys.argv[1])
manifest=json.loads((root/"manifest.json").read_text())
output=root/"work"/"outputs"
gallery=output/"showcase"/"heat_wave_annulus"
rc=int((root/"native.exitcode").read_text())
status={"schema":"v130-heat-wave-long-native-visuals-v1",
        "status":"FAIL_OR_UNQUALIFIED_LONG_HEAT_WAVE_SHOWCASE",
        "native_exit_code":rc,"profile":manifest["profile"],
        "true_FE_long_physics_gates_qualified":False,
        "whispering_gallery_mode_qualified":False,"errors":[]}
try:
    if rc: raise RuntimeError("Long native GetFEM/PETSc failed")
    report=json.loads((gallery/"provenance.json").read_text())
    with np.load(gallery/"heat_wave_FE_H_long.npz",allow_pickle=False) as arr:
        t=arr["t"]; H=arr["H"]
    if not (report["schema"]=="scrimp-heat-wave-long-native-showcase-v1"
            and report["source"].startswith("examples/heat_wave.py")
            and report["initial_profile"]==manifest["profile"]
            and report["no_time_interpolation"] is True
            and not report["whispering_gallery_eigenmode_certified"]
            and not report["full_field_energy_localization_certified"]
            and H.shape==t.shape and len(t)>100
            and abs(t[0])<1e-12 and abs(t[-1]-15.)<1e-9
            and np.all(np.diff(t)>0.) and np.all(np.isfinite(H))
            and np.all(H>=0.) and H[0]>0.
            and abs(H[0]-report["initial_H"])<1e-11
            and abs(H[-1]-report["final_H"])<1e-11):
        raise RuntimeError("Invalid real long-time FE history or gallery caveats")
    if float(np.max(np.diff(H)))>1e-7*max(1.,float(H[0])):
        raise RuntimeError("Long physical Hamiltonian has unaccounted energy growth")
    files=["heat_wave_H_long_loglog.png",
           "heat_wave_H_long_loglog.svg",
           "wave_annulus_long.gif",
           "wave_outer_shell_kinetic_proxy.png"]
    for name in files:
        f=gallery/name
        if not f.is_file() or f.stat().st_size<250:
            raise RuntimeError("Missing in-example visual "+name)
        signature=f.read_bytes()[:6]
        if name.endswith(".png") and not signature.startswith(b"\x89PNG"):
            raise RuntimeError("Invalid native PNG "+name)
        if name.endswith(".gif") and signature not in (b"GIF89a",b"GIF87a"):
            raise RuntimeError("Invalid native GIF")
        if name.endswith(".svg"):
            ET.parse(f)
    with Image.open(gallery/"wave_annulus_long.gif") as im:
        if im.n_frames!=report["gif_frames_real"]:
            raise RuntimeError("The gif omitted original GetFEM FE states")
    for name in ("p","q","T"):
        pvd=output/"pv"/name/(name+".pvd")
        frames=ET.parse(pvd).findall(".//DataSet")
        if len(frames)!=len(t):
            raise RuntimeError("Missing long native VTK frames for "+name)
        for frame in frames:
            rel=Path(frame["file"])
            if rel.name!=str(rel) or not (pvd.parent/rel).is_file():
                raise RuntimeError("Missing native GetFEM VTU for "+name)
    status["status"]="PASS_LONG_NATIVE_VISUALS_NOT_WG_MODE_PROOF"
    status["times_saved"]=int(len(t))
    status["gif_frames"]=int(report["gif_frames_real"])
    status["max_H_increment"]=float(np.max(np.diff(H)))
    status["H_initial"]=float(H[0])
    status["H_final"]=float(H[-1])
except Exception as exc:
    status["errors"].append(f"{type(exc).__name__}: {exc}")
(root/"summary.json").write_text(json.dumps(status,indent=2)+"\n")
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,compresslevel=4,
                     allowZip64=True) as dest:
    for f in sorted(root.rglob("*")):
        if f.is_file():dest.write(f,arcname=str(f.relative_to(root)))
print(json.dumps({**status,"archive":str(archive),
                  "evidence":str(root)},indent=2),flush=True)
raise SystemExit(0 if status["status"].startswith("PASS_") else 1)
PY
