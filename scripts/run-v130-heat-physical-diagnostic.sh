#!/usr/bin/env bash
# SCRIMP v1.3.0: full original heat example, independent native FE/power traces.
# The original model is executed unchanged; this is a diagnostic, not physics PASS.
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
  echo "Unexpected native Docker scientific revision: $revision" >&2
  exit 2
}
python3 -m py_compile examples/heat.py sandbox/probe_v130_heat_physical.py
bash -n "$0"
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-heat-physical-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
echo "SCRIMP v1.3.0 original 2D heat — full FE + original GetFEM physical power diagnostic"
echo "Canonical 500 BDF2 steps, unchanged h=0.1, dt=0.01, t_f=5; MPI=1/2"
echo "EVIDENCE: $out"
printf '{"schema":"scrimp-v130-heat-original-physical-diagnostic-v1","source":"%s","image_revision":"%s","model":"examples/heat.py","MPI":[1,2],"time_scheme":"original BDF2","physical_PASS_not_claimed":true}\n' \
  "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
all_rc=0
for n in 1 2; do
  mkdir -p "$out/native_mpi$n" "$out/work_mpi$n"
  echo "RUN native MPI=$n"
  case_rc=0
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
    --volume "$root/examples/heat.py:/workspace/scrimp/examples/heat.py:ro" \
    --volume "$root/sandbox/probe_v130_heat_physical.py:/workspace/scrimp/sandbox/probe_v130_heat_physical.py:ro" \
    --volume "$out/native_mpi$n:/artifacts:rw" \
    --volume "$out/work_mpi$n:/scrimp-work:rw" \
    --workdir /scrimp-work \
    "$image" timeout -k 10s "${SCRIMP_V130_HEAT_TIMEOUT:-1800s}" \
    mpiexec --oversubscribe --bind-to none -n "$n" \
    python -u /workspace/scrimp/sandbox/probe_v130_heat_physical.py \
      --out /artifacts \
    > "$out/mpi$n.host.log" 2>&1 || case_rc=$?
  printf "%s\n" "$case_rc" > "$out/mpi$n.exitcode"
  if [[ "$case_rc" -ne 0 ]]; then
    all_rc=1
    echo "UNQUALIFIED native original heat MPI=$n exit=$case_rc" >&2
    tail -n 65 "$out/mpi$n.host.log" >&2
  else
    echo "NATIVE COMPLETE MPI=$n — original physical traces stored"
  fi
done

compare_rc=0
python3 - "$out" <<'PY' || compare_rc=$?
import json
import sys
from pathlib import Path
import numpy as np

root=Path(sys.argv[1])
problems=[]
reports={}
datasets={}
for n in (1, 2):
    if (root/f"mpi{n}.exitcode").read_text().strip()!="0":
        problems.append(f"Native MPI{n} command did not finish")
        continue
    wanted=(0,) if n==1 else (0,1)
    for rank in wanted:
        key=(n,rank)
        p=root/f"native_mpi{n}/rank{rank}.json"
        q=root/f"native_mpi{n}/rank{rank}.npz"
        if not p.is_file() or not q.is_file():
            problems.append(f"Original native trace missing: MPI{n}, rank {rank}")
            continue
        reports[key]=json.loads(p.read_text())
        with np.load(q) as f:
            datasets[key]={k:np.array(f[k]) for k in f.files}
    fig=root/f"work_mpi{n}/outputs/png/Hamiltonian_Heat_2D.png"
    if not fig.is_file() or fig.stat().st_size<1000:
        problems.append(f"Original Hamiltonian PNG missing for MPI{n}")

comparisons={}
if (1,0) in datasets and (2,0) in datasets and (2,1) in datasets:
    ref=datasets[1,0]
    for key in ((2,0),(2,1)):
        x=datasets[key]
        comparisons[f"mpi1_vs_mpi{key[0]}_rank{key[1]}"]={
            "max_complete_FE_state_absolute_difference":float(np.max(
                np.abs(ref["states"]-x["states"]))),
            "max_H_absolute_difference":float(np.max(
                np.abs(ref["H"]-x["H"]))),
            "max_Fourier_dissipation_absolute_difference":float(np.max(
                np.abs(ref["dissipation"]-x["dissipation"]))),
            "max_outward_boundary_flux_absolute_difference":float(max(
                np.max(np.abs(ref[f"boundary_{r}"]-x[f"boundary_{r}"]))
                for r in (10,11,12,13)
            )),
            "max_correct_port_power_absolute_difference":float(
                np.max(np.abs(ref["control_power"]-x["control_power"]))),
            "max_BDF2_quadratic_defect_MPI_absolute_difference":float(
                np.max(np.abs(
                    ref["discrete_bdf2_power_defect"][2:]
                    - x["discrete_bdf2_power_defect"][2:]))),
        }
    for key, data in datasets.items():
        if (data["states"].shape!=(501,5830)
                or data["H"].shape!=(501,)
                or not np.all(np.isfinite(data["states"]))):
            problems.append(f"Invalid complete native data shape/finite: {key}")
report={
    "schema":"scrimp-v130-heat-physical-host-comparison-v1",
    "status":("NATIVE_COMPLETE_PHYSICS_REVIEW_PENDING" if not problems
              else "FAIL_OR_INCOMPLETE_NATIVE_PHYSICS_DIAGNOSTIC"),
    "errors":problems,
    "not_an_exact_BDF2_energy_identity":True,
    "physical_energy_acceptance_pending_independent_review":True,
    "MPI_full_FE_and_power_comparisons":comparisons,
    "native_rank_reports":{
        f"mpi{k[0]}_rank{k[1]}":v for k,v in reports.items()
    },
}
(root/"summary.json").write_text(json.dumps(report,indent=2)+"\n")
print(json.dumps({**report,"native_rank_reports":{
    k:{
        "H_initial":v["H_initial"],"H_final":v["H_final"],
        "max_BDF_trapezoid_energy_defect":v[
            "max_BDF_trapezoid_energy_defect"
        ],
        "relative_BDF_trapezoid_defect":v[
            "relative_BDF_trapezoid_defect"
        ],
        "max_H_independent_difference":v["max_H_independent_difference"],
        "max_constitutive_dissipation_difference":v[
            "max_constitutive_dissipation_difference"
        ],
        "max_correct_port_trapezoid_energy_defect":v[
            "max_correct_port_trapezoid_energy_defect"
        ],
        "max_exact_BDF2_quadratic_power_defect":v[
            "max_exact_BDF2_quadratic_power_defect"
        ],
        "BDF1_initial_step_power_defect":v[
            "BDF1_initial_step_power_defect"
        ],
    } for k,v in report["native_rank_reports"].items()
}},indent=2),flush=True)
raise SystemExit(bool(problems))
PY
if [[ "$compare_rc" -ne 0 ]]; then all_rc=1; fi

python3 - "$out" <<'PY'
import json
import sys
import zipfile
from pathlib import Path
root=Path(sys.argv[1]); archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                     compresslevel=4,allowZip64=True) as zf:
    for f in sorted(root.rglob("*")):
        if f.is_file():
            zf.write(f,arcname=str(f.relative_to(root)))
print(json.dumps({"archive":str(archive),
                  "evidence":str(root),
                  "summary":str(root/"summary.json")},indent=2))
PY
exit "$all_rc"
