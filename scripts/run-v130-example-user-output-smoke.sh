#!/usr/bin/env bash
# End-user contract: each canonical v1.3 example must DIRECTLY create
# its own real visual outputs under work/outputs/<system_name>/<kind>.
# This bounded test does NOT replace already qualified long MPI physics.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
    --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
if [[ "$revision" != 1d6d0c90115213298047b28237a7f38087563a35 ]]; then
  echo "Unqualified GetFEM/PETSc rc3 image revision $revision" >&2
  exit 2
fi
python3 tests/test_v130_user_output_power_hamiltonian.py
python3 tests/test_v130_examples_output_layout_source.py
python3 tests/test_v130_nonlinear_pendulum_source.py
python3 scripts/check-v130-published-immutable.py
python3 -m py_compile \
  scrimp/dphs.py scrimp/port.py scrimp/hamiltonian.py \
  scrimp/utils/example_output.py scrimp/utils/power_sampling.py \
  examples/heat.py examples/heat_wave.py \
  examples/mass_spring_damper.py examples/nonlinear_pendulum.py \
  examples/shallow_water.py examples/wave.py \
  examples/wave_dissipative.py examples/wave_coenergy.py
bash -n "$0"

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-user-examples-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
echo "SCRIMP v1.3.0 original example outputs — real native GetFEM, NOT a source-only mock"
echo "All shared images and ParaView files in runs/<case>/outputs/<system>/{png,svg,gif,mesh,log,npz,json}"
echo "Docker image $image, revision $revision"
echo "Evidence directory: $out"

# Calls the original public Python example entry point; every case gets
# a separate actual user working_directory. MPI processes share ONE
# output root, with ONE validated shared ParaView writer.
run_case() {
    local label="$1" system="$2" mpi="$3" module="$4"
    shift 4
    local wd="$out/runs/$label"
    local reportdir="$wd/outputs/$system/log"
    mkdir -p "$reportdir"
    echo "RUN canonical $module, case=$label, MPI=$mpi, user system=$system"
    local rc=0
    "${docker_cmd[@]}" run --rm \
        --user "$(id -u):$(id -g)" \
        --env HOME=/tmp --env MPLBACKEND=Agg \
        --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
        --env PYTHONPATH=/workspace/scrimp \
        --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
        --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
        --env SCRIMP_MATRIX_MEMORY_MODE=compact \
        --env SCRIMP_MUMPS_RESET_POLICY=on_change \
        --env SCRIMP_GETFEM_RHS_MODE=auto \
        --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
        --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
        --volume "$root/examples:/workspace/scrimp/examples:ro" \
        --volume "$wd:/scrimp-work:rw" \
        --workdir /scrimp-work "$image" \
        timeout -k 10s "${SCRIMP_V130_EXAMPLE_SMOKE_TIMEOUT:-900s}" \
        mpiexec --oversubscribe --bind-to none -n "$mpi" \
        python -u -m "$module" "$@" \
        >"$reportdir/run.log" 2>&1 || rc=$?
    printf '%s\n' "$rc" > "$reportdir/exitcode.txt"
    if [[ "$rc" != 0 ]]; then
        echo "FAIL native user example: $label (exit=$rc)" >&2
        tail -n 90 "$reportdir/run.log" >&2
        echo "Evidence retained: $out" >&2
        exit "$rc"
    fi
    python3 - "$wd" "$system" "$label" <<'PY'
from pathlib import Path
import json,sys
from xml.etree import ElementTree as ET
work=Path(sys.argv[1]); system=sys.argv[2]; label=sys.argv[3]
root=work/"outputs"/system
png=list((root/"png").glob("*.png"))
svg=list((root/"svg").glob("*.svg"))
data=list((root/"npz").glob("*.npz"))
pv=list((root/"mesh").glob("*/*.pvd"))
gif=list((root/"gif").glob("*.gif"))
if not png or not svg or not data:
    raise RuntimeError(f"{label}: missing canonical PNG/SVG/NPZ in {root}")
for p in [*png,*svg,*data,*gif,*pv]:
    if not p.is_file() or p.stat().st_size<100:
        raise RuntimeError(f"{label}: invalid native output {p}")
for p in pv:
    sets=ET.parse(p).getroot().findall(".//DataSet")
    if len(sets)<2:
        raise RuntimeError(f"{label}: missing true GetFEM VTK time sequence {p}")
    times=[float(v.attrib["timestep"]) for v in sets]
    if any(b<=a for a,b in zip(times,times[1:])):
        raise RuntimeError(f"{label}: invalid true GetFEM PVD clock {p}")
    for obj in sets:
        f=p.parent/obj.attrib["file"]
        if not f.is_file():
            raise RuntimeError(f"{label}: missing true GetFEM VTU {f}")
manifest={
    "case":label,"status":"PASS_CANONICAL_REAL_GETFEM_USER_OUTPUTS",
    "working_directory":str(work),"system_name":system,
    "output_root":str(root),"png":len(png),"svg":len(svg),
    "npz":len(data),"gif":len(gif),"pvd":len(pv),
    "readable_true_VTK_clocks":True if pv else None,
}
print(json.dumps(manifest,indent=2),flush=True)
(root/"json").mkdir(exist_ok=True)
(root/"json"/"user_output_smoke.json").write_text(
    json.dumps(manifest,indent=2)+"\n",encoding="utf-8"
)
PY
    echo "PASS $label: $wd/outputs/$system/"
}

# The bounded runs are only for CURRENT OUTPUT and public API regression,
# not a repeat of all historical full-horizon scientific qualifications.
run_case mass_spring mass_spring_damper 1 examples.mass_spring_damper
run_case pendulum nonlinear_pendulum_free 1 examples.nonlinear_pendulum
run_case heat heat 1 examples.heat --final-time .06 --time-step .01
run_case heat_wave heat_wave 1 examples.heat_wave \
    --final-time .02 --time-step .002 --save-step .002
run_case shallow_water shallow_water 1 examples.shallow_water \
    --final-time .02 --time-step .0001 --save-step .002
run_case wave_2d wave_2d 1 examples.wave
run_case wave_dissipative wave_dissipative 1 examples.wave_dissipative
run_case wave_coenergy wave_coenergy 1 examples.wave_coenergy
# Specific MPI shared-writer regression for the new output path;
# physics MPI1/2 long-horizon is already separately accepted.
run_case heat_wave_mpi2 heat_wave 2 examples.heat_wave \
    --final-time .02 --time-step .002 --save-step .002

python3 - "$out" "$image" "$revision" <<'PY'
import json,subprocess,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1])
reports=sorted(root.glob("runs/*/outputs/*/json/user_output_smoke.json"))
if len(reports)!=9:
    raise RuntimeError(f"Missing native canonical output checks: {len(reports)} != 9")
result={"status":"PASS_V130_ALL_NINE_REAL_USER_EXAMPLE_OUTPUT_CASES",
 "schema":"scrimp-v130-all-user-canonical-output-layout-v1",
 "source":subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip(),
 "image":sys.argv[2],"image_revision":sys.argv[3],
 "cases":[json.loads(p.read_text()) for p in reports],
 "no_scientific_order_or_broad_MPI_convergence_claim":True}
(root/"summary.json").write_text(json.dumps(result,indent=2)+"\n")
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",zipfile.ZIP_DEFLATED,compresslevel=3) as z:
    for p in sorted(root.rglob("*")):
        if p.is_file():
            z.write(p,str(p.relative_to(root)))
with zipfile.ZipFile(archive) as z:
    bad=z.testzip()
    if bad:
        raise RuntimeError(f"Corrupted native output archive: {bad}")
print(json.dumps({"status":result["status"],"evidence":str(root),
                  "archive":str(archive),"ZIP_CRC_pass":True},
                 indent=2),flush=True)
PY
