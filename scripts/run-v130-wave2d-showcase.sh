#!/usr/bin/env bash
# SCRIMP v1.3.0: real 2D wave showcase (2A, 2B or 2C), one MUMPS rank.
# Longer native GetFEM/PETSc trajectory for user-facing visual quality.
# This does NOT replace the short MUMPS MPI1/2 scientific qualification.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 1 && ( $1 == 2a || $1 == 2b || $1 == 2c ) ]] || {
    echo "Usage: $0 {2a|2b|2c}" >&2; exit 2;
}
case "$1" in
    2a) module=wave;;
    2b) module=wave_dissipative;;
    2c) module=wave_coenergy;;
esac
docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
  --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
    echo "Wrong pinned native GetFEM/PETSc image $revision" >&2; exit 2;
}
h="${SCRIMP_V130_SHOWCASE_H:-0.1}"
tf="${SCRIMP_V130_SHOWCASE_TF:-1.0}"
dt="${SCRIMP_V130_SHOWCASE_DT:-0.01}"
python3 -m py_compile "examples/$module.py" scripts/render-v130-wave2d-native.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-wave2d-$1-showcase-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
mkdir -p "$out/mumps_mpi1" "$out/work"
printf '{"schema":"scrimp-v130-native-wave2d-long-visual-showcase-v1","source":"%s","image":"%s","image_revision":"%s","canonical":"examples/%s.py","scenario":"%s","mesh_size":"%s","final_time":"%s","time_step":"%s","MPI":1,"solver":"mumps","not_an_MPI_parity_or_convergence_qualification":true}\n' \
  "$(git rev-parse HEAD)" "$image" "$revision" "$module" "$1" "$h" "$tf" "$dt" >"$out/manifest.json"
echo "Real SCRIMP 2D $1 showcase: h=$h, final_time=$tf, dt=$dt"
echo "Will export native GetFEM FE, VTK, CSV, PNG, SVG and animated GIF."
echo "EVIDENCE: $out"
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
    --volume "$root/examples/$module.py:/workspace/scrimp/examples/$module.py:ro" \
    --volume "$out/mumps_mpi1:/artifacts:rw" \
    --volume "$out/work:/scrimp-work:rw" \
    --workdir /scrimp-work \
    "$image" timeout -k 5s "${SCRIMP_V130_SHOWCASE_TIMEOUT:-1200s}" \
    mpiexec --oversubscribe --bind-to none -n 1 \
    python -u -m "examples.$module" \
       --solver mumps --mesh-size "$h" --final-time "$tf" --time-step "$dt" \
       --check --export --output /artifacts \
    >"$out/mumps_mpi1.host.log" 2>&1 || rc=$?
echo "$rc" >"$out/mumps_mpi1.exitcode"
if [[ "$rc" == 0 ]]; then
    python3 scripts/render-v130-wave2d-native.py \
      --input "$out" --case mumps_mpi1 --output "$out/visuals" || rc=$?
else
    echo "Unqualified native long showcase, exit=$rc" >&2
    tail -n 55 "$out/mumps_mpi1.host.log" >&2
fi
python3 - "$out" "$rc" <<'PY'
import json,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1]);rc=int(sys.argv[2])
summary={"schema":"v130-wave2d-showcase-v1",
         "status":("PASS_NATIVE_LONG_VISUAL_SHOWCASE"
                   if rc==0 else "FAIL_OR_UNQUALIFIED_SHOWCASE"),
         "native_or_visual_exit_code":rc,
         "not_MPI_parity_or_time_mesh_convergence":True}
(root/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as out:
    for file in sorted(root.rglob("*")):
        if file.is_file():
            out.write(file,arcname=str(file.relative_to(root)))
print(json.dumps({**summary,"archive":str(archive),
                  "outputs":str(root)},indent=2))
raise SystemExit(0 if rc==0 else 1)
PY
