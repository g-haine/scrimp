#!/usr/bin/env bash
# SCRIMP 1.3.0 — native 4-step acceptance pilot for improved canonical heat.
# The baseline BDF2/pseudo model remains independently reproduced by
# scripts/run-v130-heat-physical-diagnostic.sh.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ "$#" == 0 ]] || { echo "Usage: $0" >&2; exit 2; }

docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
    --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
    echo "Unqualified GetFEM/PETSc Docker revision: $revision" >&2
    exit 2
}
python3 -m py_compile examples/heat.py
python3 tests/test_v130_heat_cn_source.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-heat-CN-pilot-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-heat-CN-algebraic-4-step-pilot-v1","source":"%s","image":"%s","image_revision":"%s","original_FEM_and_controls_unchanged":true,"time_stepper":"cn","initializer":"algebraic","MPI":[1,2],"CN_steps":4,"full_run_and_visuals_not_yet_qualified":true}\n' \
    "$(git rev-parse HEAD)" "$image" "$revision" > "$out/manifest.json"

echo "SCRIMP 2-D heat CN + exact algebraic DAE projection — native four-step pilot"
echo "Original physical mesh h=.1, dt=.01, t_f=.04; MUMPS, COO; MPI1/2"
echo "No visual showcase yet — mathematical acceptance first"
echo "EVIDENCE: $out"
overall=0
for n in 1 2; do
    mkdir -p "$out/native_mpi$n" "$out/work_mpi$n"
    echo "RUN CN/algebraic MPI=$n"
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
        --volume "$root/examples/heat.py:/workspace/scrimp/examples/heat.py:ro" \
        --volume "$out/native_mpi$n:/artifacts:rw" \
        --volume "$out/work_mpi$n:/scrimp-work:rw" \
        --workdir /scrimp-work \
        "$image" timeout -k 10s "${SCRIMP_V130_HEAT_CN_TIMEOUT:-900s}" \
        mpiexec --oversubscribe --bind-to none -n "$n" \
        python -u -m examples.heat \
          --time-stepper cn --initializer algebraic \
          --mesh-size 0.1 --time-step 0.01 --final-time 0.04 \
          --solver mumps --no-showcase --check --output /artifacts \
        > "$out/mpi$n.host.log" 2>&1 || rc=$?
    printf '%s\n' "$rc" > "$out/mpi$n.exitcode"
    if [[ "$rc" -ne 0 ]]; then
        overall=1
        echo "UNQUALIFIED improved heat CN MPI=$n exit=$rc" >&2
        tail -n 80 "$out/mpi$n.host.log" >&2
    else
        echo "NATIVE COMPLETE improved heat CN MPI=$n"
    fi
done

compare_rc=0
python3 - "$out" <<'PY' || compare_rc=$?
import json
import sys
from pathlib import Path
import numpy as np

root=Path(sys.argv[1]); errors=[]; traces={}; reports={}
for mpi in (1,2):
    if (root/f"mpi{mpi}.exitcode").read_text().strip()!="0":
        errors.append(f"Native MPI{mpi} did not complete successfully")
        continue
    fig=root/f"work_mpi{mpi}/outputs/png/Hamiltonian_Heat_2D.png"
    if not fig.is_file() or fig.stat().st_size<1000:
        errors.append(f"Original user Hamiltonian figure missing MPI{mpi}")
    for rank in range(mpi):
        base=root/f"native_mpi{mpi}/rank{rank}"
        p=Path(str(base)+".json")
        q=Path(str(base)+".npz")
        if not p.is_file() or not q.is_file():
            errors.append(f"Missing original accepted FE data MPI{mpi}, rank{rank}")
            continue
        key=(mpi,rank)
        reports[key]=json.loads(p.read_text())
        with np.load(q) as f:
            traces[key]={k:np.asarray(f[k]) for k in f.files}
        meta=reports[key]
        init=meta.get("initialization_report") or {}
        gates=(meta.get("physics") or {}).get("gates") or {}
        if (meta.get("status")!="PASS_CN_HEAT_NATIVE_PHYSICS"
                or meta.get("time_stepper")!="cn"
                or meta.get("initializer")!="algebraic"
                or meta.get("steps")!=4
                or meta.get("DOFs")!=5830
                or init.get("status")!="PASS"
                or init.get("strategy")!="algebraic"
                or init.get("differential_initial_exact") is not True
                or init.get("pseudo_ts_executed") is not False
                or not gates or not all(gates.values())):
            errors.append(f"Original native FE/DAE/CN physics gate failed {key}")
        data=traces[key]
        if (data["states"].shape!=(5,5830)
                or data["H"].shape!=(5,)
                or data["time"].shape!=(5,)
                or not np.all(np.isfinite(data["states"]))
                or np.max(np.abs(data["time"]-np.linspace(0.,.04,5)))>1e-10):
            errors.append(f"Unexpected accepted genuine CN trajectory {key}")

parity={}
if len(traces)==3:
    base=traces[(1,0)]
    for key in ((2,0),(2,1)):
        data=traces[key]
        fe=float(np.max(np.abs(base["states"]-data["states"])))
        h=float(np.max(np.abs(base["H"]-data["H"])))
        parity[f"MPI1_to_MPI{key[0]}_rank{key[1]}"]={
            "max_FE_coefficient_abs_diff":fe,
            "max_H_abs_diff":h,
        }
        if fe>1e-8 or h>1e-9:
            errors.append(f"Native MPI FE or H numerical mismatch: {key}")
summary={
    "schema":"scrimp-v130-heat-CN-algebraic-4-step-pilot-v1",
    "status":("PASS_NATIVE_HEAT_CN_ALGEBRAIC_MPI1_MPI2_PILOT"
              if not errors else "FAIL_OR_UNQUALIFIED_HEAT_CN_PILOT"),
    "errors":errors,
    "four_step_short_pilot_only":True,
    "full_original_horizon_and_showcase_pending":True,
    "full_FE_MPI_parity":parity,
    "ranks":{
        f"mpi{n}_rank{k}":{
            "init": reports[(n,k)].get("initialization_report"),
            "physical": reports[(n,k)].get("physics"),
        } for n,k in reports
    },
}
(root/"summary.json").write_text(json.dumps(summary,indent=2,default=str)+"\n")
print(json.dumps(summary,indent=2,default=str),flush=True)
raise SystemExit(bool(errors))
PY
if [[ "$compare_rc" -ne 0 ]]; then overall=1; fi

python3 - "$out" <<'PY'
import json
import sys
import zipfile
from pathlib import Path

root=Path(sys.argv[1]); archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
    for f in sorted(root.rglob("*")):
        if f.is_file(): z.write(f,arcname=str(f.relative_to(root)))
print(json.dumps({"archive":str(archive),"evidence":str(root),
                  "summary":str(root/"summary.json")},indent=2))
PY
exit "$overall"
