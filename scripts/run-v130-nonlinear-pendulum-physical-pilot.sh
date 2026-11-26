#!/usr/bin/env bash
# SCRIMP v1.3 original nonlinear lumped pendulum: SOURCE-STAGED native pilot.
# Full sine restoring torque, MUMPS MPI1/2, unforced and forced cases.
# Keep original GetFEM bricks and BDF2 integrator; exact IC by default.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }
docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
  --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == 1d6d0c90115213298047b28237a7f38087563a35 ]] || {
  echo "Unqualified Docker GetFEM/PETSc revision: $revision" >&2; exit 2;
}

python3 -m py_compile \
  examples/nonlinear_pendulum.py \
  sandbox/probe_v130_nonlinear_pendulum_physics.py \
  tests/test_v130_nonlinear_pendulum_source.py
python3 tests/test_v130_nonlinear_pendulum_source.py
python3 tests/test_v130_finite_dimensional_source.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-nonlinear-pendulum-physical-pilot-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
python3 - "$out" "$image" "$revision" <<'PY'
import json,subprocess,sys
from pathlib import Path
out=Path(sys.argv[1])
(out/"manifest.json").write_text(json.dumps({
    "schema":"scrimp-v130-original-nonlinear-pendulum-short-pilot-v1",
    "source":subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip(),
    "image":sys.argv[2],
    "image_revision":sys.argv[3],
    "original_example":"examples/nonlinear_pendulum.py",
    "source_staged_native_pending":True,
    "original_full_sine_model_unchanged":True,
    "original_BDF2_unchanged":True,
    "exact_coenergy_initial_without_pseudo":True,
    "original_angle0":1.2,
    "original_momentum0":0.,
    "cases":[
        {"name":"free","torque":0.,"MPI":[1,2]},
        {"name":"forced","torque":.1,"MPI":[1,2]},
    ],
    "t_f":.08,"dt":.01,
    "nonquadratic_power_identity":
        "e_theta D2theta + e_p D2p + d e_p² - tau e_p = 0",
    "no_false_quadratic_H_BDF2_claim":True,
},indent=2)+"\n")
PY
echo "SCRIMP v1.3 nonlinear pendulum ORIGINAL full-sine pH/DAE short physical pilot"
echo "Exact (theta0,p0,e_theta0,e_p0)=(1.2,0,9.81 sin(1.2),0)"
echo "Nonlinear energy/effort, BDF2 dt=.01, t_f=.08, MUMPS MPI1+MPI2"
echo "Two genuine original force cases: unforced and torque amplitude 0.1"
echo "True nonquadratic BDF2 coenergy/ODE power, independent full-sine RK4"
echo "No code inside examples/published modified; no figures or full-horizon claim"
echo "EVIDENCE: $out"

overall=0
for case in "free:0" "forced:0.1"; do
    IFS=: read -r name torque <<< "$case"
    for mpi in 1 2; do
        label="${name}_mpi${mpi}"
        mkdir -p "$out/native_$label" "$out/work_$label"
        echo "RUN original nonlinear pendulum $label torque=$torque"
        rc=0
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
            --volume "$root/examples/nonlinear_pendulum.py:/workspace/scrimp/examples/nonlinear_pendulum.py:ro" \
            --volume "$root/sandbox/probe_v130_nonlinear_pendulum_physics.py:/probe-pend.py:ro" \
            --volume "$out/native_$label:/artifacts:rw" \
            --volume "$out/work_$label:/scrimp-work:rw" \
            --workdir /scrimp-work \
            "$image" timeout -k 10s "${SCRIMP_V130_PENDULUM_PILOT_TIMEOUT:-420s}" \
            mpiexec --oversubscribe --bind-to none -n "$mpi" \
            python -u /probe-pend.py \
                --final-time .08 --dt .01 --angle 1.2 --torque "$torque" \
                --solver mumps \
                --out /artifacts \
            > "$out/$label.host.log" 2>&1 || rc=$?
        printf '%s\n' "$rc" > "$out/$label.exitcode"
        if [[ "$rc" != 0 ]]; then
            echo "UNQUALIFIED original nonlinear $label exit=$rc" >&2
            tail -n 90 "$out/$label.host.log" >&2
            overall=1
        else
            echo "NATIVE FINISHED: $label"
        fi
    done
done

compare_rc=0
python3 - "$out" <<'PY' || compare_rc=$?
import json,sys,zipfile
from pathlib import Path
import numpy as np
root=Path(sys.argv[1])
errors=[]
reports={}
traces={}
parity={}
for case,torque in (("free",0.),("forced",.1)):
    for mpi in (1,2):
        label=f"{case}_mpi{mpi}"
        code=root/f"{label}.exitcode"
        if not code.exists() or code.read_text().strip()!="0":
            errors.append(f"Native original nonlinear run failed: {label}")
            continue
        for rank in range(mpi):
            file=root/f"native_{label}/rank{rank}"
            try:
                info=json.loads(Path(str(file)+".json").read_text())
                with np.load(Path(str(file)+".npz"),allow_pickle=False) as z:
                    arr={k:np.asarray(z[k]) for k in z.files}
                key=(case,mpi,rank)
                reports[key]=info
                traces[key]=arr
                if (info["status"]!="PASS_ORIGINAL_NONLINEAR_PENDULUM_NATIVE_PILOT"
                        or info["steps"]!=8 or info["saved_states"]!=9
                        or info["DAE_DOFs"]!=4 or info["initialization"]!="exact_constitutive"
                        or info["pseudo_TS_executed"] is not False
                        or info["original_fully_nonlinear_sine"] is not True
                        or abs(info["torque_amplitude"]-torque)>1e-14
                        or not all(info["gates"].values())
                        or arr["t"].shape!=(9,)
                        or arr["z"].shape!=(9,4)
                        or not np.all(np.isfinite(arr["z"]))
                        or np.max(np.abs(arr["t"]-np.linspace(0.,.08,9)))>1e-10):
                    errors.append(f"Original nonlinear physical gate failed: {key}")
            except Exception as exc:
                errors.append(f"Cannot read native original pendulum {label}/{rank}: {exc}")
    if all((case,mpi,rank) in traces for mpi,rank in ((1,0),(2,0),(2,1))):
        ref=traces[case,1,0]
        for rank in (0,1):
            other=traces[case,2,rank]
            discrepancies={}
            for name in ("t","z","theta","p","e_theta","e_p","H",
                         "torque","input_power","viscous_dissipation",
                         "BDF2_nonquadratic_coenergy_power_defect"):
                a,b=ref[name],other[name]
                if a.shape!=b.shape:
                    errors.append(f"Original pendulum MPI shape {case}/{rank}/{name}")
                    continue
                delta=float(np.max(np.abs(a-b)))
                discrepancies[name]=delta
                if not np.isfinite(delta) or delta>1e-8:
                    errors.append(f"Original pendulum MPI mismatch {case}/{rank}/{name}: {delta}")
            parity[f"{case}_mpi1_to_mpi2_rank{rank}"]=discrepancies
        if not np.array_equal(traces[case,2,0]["z"],traces[case,2,1]["z"]):
            errors.append(f"Original nonlinear MPI2 replicas differ: {case}")

summary={
    "schema":"scrimp-v130-original-nonlinear-pendulum-native-short-v1",
    "status":("PASS_ORIGINAL_PENDULUM_FULL_SINE_BDF2_FORCED_UNFORCED_MPI1_MPI2"
              if not errors else "FAIL_OR_UNQUALIFIED_PENDULUM_SHORT_PILOT"),
    "errors":errors,
    "full_sine_nonlinear_constitutive":True,
    "no_false_quadratic_BDF2_H_identity":True,
    "fully_preserved_physical_initial_data":not errors,
    "short_pilot_only":True,
    "full_horizon_and_canonical_visuals_pending":True,
    "parity_all_four_original_lumped_DAE_variables":parity,
    "native_reports":{f"{case}_mpi{mpi}_rank{rank}":meta
                      for (case,mpi,rank),meta in reports.items()},
}
(root/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
    for file in sorted(root.rglob("*")):
        if file.is_file():
            z.write(file,arcname=str(file.relative_to(root)))
with zipfile.ZipFile(archive) as z:
    if z.testzip() is not None:
        raise RuntimeError("Nonlinear original GetFEM evidence ZIP CRC failed")
print(json.dumps({**summary,"archive":str(archive),
                  "evidence":str(root),"ZIP_CRC_pass":True},indent=2))
raise SystemExit(bool(errors))
PY
if [[ "$compare_rc" != 0 ]]; then overall=1; fi
echo "ARCHIVE: $out.zip"
if [[ "$overall" == 0 ]]; then
    echo "PASS_ORIGINAL_PENDULUM_FULL_SINE_BDF2_FORCED_UNFORCED_MPI1_MPI2"
else
    echo "FAIL_OR_UNQUALIFIED_PENDULUM_SHORT_PILOT" >&2
fi
exit "$overall"
