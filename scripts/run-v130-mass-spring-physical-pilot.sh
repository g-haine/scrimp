#!/usr/bin/env bash
# SCRIMP v1.3.0 -- original finite-dimensional linear example acceptance.
# Real GetFEM lumped port state, canonical BDF2 and exact coenergy IC, MPI1/2.
# This is a short physical pilot only; long run/figures remain separate.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }

docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == "1" ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
    --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
    echo "Unqualified GetFEM/PETSc Docker image revision: $revision" >&2
    exit 2
}

python3 -m py_compile \
    examples/mass_spring_damper.py \
    sandbox/probe_v130_mass_spring_physics.py
python3 tests/test_v130_finite_dimensional_source.py
python3 scripts/check-v130-published-immutable.py

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-mass-spring-physical-pilot-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-mass-spring-original-short-pilot-v1","source":"%s","image":"%s","image_revision":"%s","original_example":"examples/mass_spring_damper.py","original_weak_model_unchanged":true,"initialization":"exact_constitutive_no_pseudo","MPI":[1,2],"final_time":0.08,"time_step":0.01,"BDF2_startup_stage_not_assumed_to_be_BDF1":true,"full_horizon_and_showcase_pending":true}\n' \
    "$(git rev-parse HEAD)" "$image" "$revision" > "$out/manifest.json"
echo "SCRIMP 1.3.0 original MASS-SPRING-DAMPER: linear lumped pH/DAE physical pilot"
echo "Original ODE m=1.5, k=5, c=0.25, f(t)=0.2 sin(2 pi t), initial (q,p)=(.4,0)"
echo "Original GetFEM bricks, BDF2 and exact algebraic coenergies at t=0, dt=.01, t_f=.08, MPI1/2"
echo "Independent RK4 ODE; true co-energies; original H; quadratic BDF2 power identity k>=2"
echo "No changes to physics, BDF2 or GetFEM bricks; exact initial data, no pseudo-time drift"
echo "EVIDENCE: $out"

overall=0
for mpi in 1 2; do
    mkdir -p "$out/native_mpi$mpi" "$out/work_mpi$mpi"
    echo "RUN original lumped ODE / MUMPS MPI=$mpi"
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
        --volume "$root/examples/mass_spring_damper.py:/workspace/scrimp/examples/mass_spring_damper.py:ro" \
        --volume "$root/sandbox/probe_v130_mass_spring_physics.py:/probe-msd.py:ro" \
        --volume "$out/native_mpi$mpi:/artifacts:rw" \
        --volume "$out/work_mpi$mpi:/scrimp-work:rw" \
        --workdir /scrimp-work \
        "$image" timeout -k 10s "${SCRIMP_V130_MSD_PILOT_TIMEOUT:-600s}" \
        mpiexec --oversubscribe --bind-to none -n "$mpi" \
        python -u /probe-msd.py \
          --final-time 0.08 --time-step 0.01 \
          --solver mumps --out /artifacts \
        > "$out/mpi$mpi.host.log" 2>&1 || rc=$?
    echo "$rc" > "$out/mpi$mpi.exitcode"
    if [[ "$rc" -ne 0 ]]; then
        overall=1
        echo "UNQUALIFIED original mass-spring MPI=$mpi exit=$rc" >&2
        tail -n 85 "$out/mpi$mpi.host.log" >&2
    else
        echo "NATIVE COMPLETE original mass-spring MPI=$mpi"
    fi
done

compare_rc=0
python3 - "$out" <<'PY' || compare_rc=$?
import json
import sys
import zipfile
from pathlib import Path
import numpy as np

root=Path(sys.argv[1])
errors=[]
reports={}
arrays={}
for mpi in (1,2):
    exitfile=root/f"mpi{mpi}.exitcode"
    if not exitfile.exists() or exitfile.read_text().strip()!="0":
        errors.append(f"Native MPI{mpi} original model did not run")
        continue
    for rank in range(mpi):
        key=(mpi,rank)
        meta=root/f"native_mpi{mpi}/rank{rank}.json"
        data=root/f"native_mpi{mpi}/rank{rank}.npz"
        if not meta.is_file() or not data.is_file():
            errors.append(f"Missing original GetFEM trace/report: MPI{mpi} rank{rank}")
            continue
        try:
            report=json.loads(meta.read_text())
            with np.load(data,allow_pickle=False) as f:
                trace={k:np.asarray(f[k]) for k in f.files}
            reports[key]=report
            arrays[key]=trace
            t=trace["t"]
            z=trace["z"]
            if (report["status"]!="PASS_ORIGINAL_LUMPED_ODE_NATIVE_PILOT"
                    or report["time_stepper"]!="original_BDF2"
                    or report["initializer"]!="exact_constitutive"
                    or report["no_pseudo_time_initialization"] is not True
                    or report["steps"]!=8
                    or report["saved_states"]!=9
                    or report["full_DAE_DOFs"]!=4
                    or not all(report["gates"].values())
                    or report["force_power_not_silently_omitted"] is not True
                    or t.shape!=(9,) or z.shape!=(9,4)
                    or not np.all(np.isfinite(z))
                    or np.max(np.abs(t-np.linspace(0.,.08,9)))>1e-10):
                errors.append(f"Original physical BDF2/ODE gates failed {key}")
        except Exception as exc:
            errors.append(f"Unreadable native GetFEM original lumped trace {key}: {exc}")
parity={}
if len(arrays)==3:
    ref=arrays[1,0]
    for key in ((2,0),(2,1)):
        other=arrays[key]
        diffs={}
        for name in ("t","z","q","p","e_q","e_p","H",
                     "external_force","external_power",
                     "Fourier_or_viscous_dissipation"):
            a,b=ref[name],other[name]
            if a.shape!=b.shape:
                errors.append(f"Full original GetFEM MPI shape mismatch {name} {key}")
                continue
            diff=float(np.max(np.abs(a-b)))
            diffs[name]=diff
            if not np.isfinite(diff) or diff>1e-8:
                errors.append(f"Full original GetFEM MPI numerical mismatch {name}: {diff}")
        parity[f"mpi1_to_mpi2_rank{key[1]}"]=diffs
    if not np.array_equal(arrays[2,0]["z"],arrays[2,1]["z"]):
        errors.append("MPI2 original full-DAE state replicas differ")
summary={
    "schema":"scrimp-v130-original-lumped-mass-spring-pH-pilot-v1",
    "status":("PASS_NATIVE_ORIGINAL_LUMPED_MASS_SPRING_MPI1_MPI2_SHORT_PILOT"
              if not errors else "FAIL_OR_UNQUALIFIED_LUMPED_PHYSICAL_PILOT"),
    "errors":errors,
    "short_pilot_only":True,
    "scientific_full_horizon_and_user_visuals_pending":True,
    "MPI1_MPI2_all_DAE_and_physical_observable_parity":parity,
    "ranks":{f"mpi{mpi}_rank{rank}":report
             for (mpi,rank),report in reports.items()},
}
(root/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2),flush=True)
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as z:
    for f in sorted(root.rglob("*")):
        if f.is_file():
            z.write(f,arcname=str(f.relative_to(root)))
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None,"ZIP CRC failure"
print(json.dumps({"archive":str(archive),"evidence":str(root),
                  "ZIP_CRC_pass":True},indent=2),flush=True)
raise SystemExit(bool(errors))
PY
if [[ "$compare_rc" -ne 0 ]]; then overall=1; fi
exit "$overall"
