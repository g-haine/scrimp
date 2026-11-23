#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# First genuine example-level PETSc LU check, not a solver benchmark.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }
docker_cmd=(docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" != 1 ]]; then
    docker_cmd=(sudo docker)
fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
    --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
    echo "Unexpected original qualified image revision: $revision" >&2
    exit 2
}

python3 -m py_compile \
    examples/mass_spring_damper.py \
    examples/nonlinear_pendulum.py \
    examples/heat_wave.py \
    examples/shallow_water.py \
    sandbox/probe_v130_lumped_example_solver.py \
    tests/test_v130_main_example_solver_policy.py
python3 tests/test_v130_main_example_solver_policy.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-example-lu-lumped-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-original-lumped-example-lu-pilot-v1","source":"%s","image_revision":"%s","status":"PENDING","physical_examples":["mass_spring_damper","nonlinear_pendulum"],"not_a_benchmark":true,"published_unchanged":true}\n' \
    "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "SCRIMP real lumped example LU qualification; originals NOT changed"
echo "MUMPS/SuperLU_DIST MPI1/2, sequential SuperLU MPI1, default MUMPS MPI1"
echo "T=0..0.04, dt=0.01, original equations and energies"
echo "EVIDENCE: $out"

for example in mass_spring_damper nonlinear_pendulum; do
    for config in default_1 mumps_1 mumps_2 superlu_1 superlu_dist_1 superlu_dist_2; do
        backend="${config%_*}"
        ranks="${config##*_}"
        case_name="${example}_${config}"
        dest="$out/$case_name"
        mkdir -p "$dest" "$out/work_$case_name"
        echo "RUN $case_name original DAE"
        if ! "${docker_cmd[@]}" run --rm \
            --user "$(id -u):$(id -g)" \
            --env HOME=/tmp --env MPLBACKEND=Agg \
            --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
            --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
            --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
            --volume "$root/examples/mass_spring_damper.py:/workspace/scrimp/examples/mass_spring_damper.py:ro" \
            --volume "$root/examples/nonlinear_pendulum.py:/workspace/scrimp/examples/nonlinear_pendulum.py:ro" \
            --volume "$root/sandbox/probe_v130_lumped_example_solver.py:/original-example-probe.py:ro" \
            --volume "$dest:/artifacts:rw" \
            --volume "$out/work_$case_name:/scrimp-work:rw" \
            --workdir /scrimp-work \
            "$image" timeout -k 5s "${SCRIMP_V130_EXAMPLE_PILOT_TIMEOUT:-180s}" \
            mpiexec --oversubscribe --bind-to none -n "$ranks" \
            python -u /original-example-probe.py \
            --example "$example" --factor-solver "$backend" \
            --out /artifacts > "$out/$case_name.host.log" 2>&1; then
            echo "FAIL original example $case_name; outputs retained $out" >&2
            tail -n 100 "$out/$case_name.host.log" >&2
            exit 1
        fi
        grep -a '^PASS ORIGINAL EXAMPLE' "$out/$case_name.host.log" || {
            echo "FAIL missing native signal $case_name" >&2
            exit 1
        }
        [[ -s "$dest/original_example_trace.npz" \
            && -s "$dest/original_example_report.json" ]] || {
            echo "FAIL original native files absent: $case_name" >&2
            exit 1
        }
    done
done

python3 - "$out" <<'PY'
import json
import sys
import zipfile
from pathlib import Path
import numpy as np

root=Path(sys.argv[1])
results={}
for example in ("mass_spring_damper","nonlinear_pendulum"):
    baseline=None
    for cfg in ("default_1","mumps_1","mumps_2",
                "superlu_1","superlu_dist_1","superlu_dist_2"):
        name=f"{example}_{cfg}"
        report=json.loads((root/name/"original_example_report.json").read_text())
        with np.load(root/name/"original_example_trace.npz",allow_pickle=False) as z:
            trace={key:z[key] for key in z.files}
        expected=cfg.rsplit("_",1)[0]
        assert report["factor_solver"]==(
            "mumps" if expected=="default" else expected
        )
        assert report["status"]=="PASS_ORIGINAL_LUMPED_EXAMPLE_NATIVE"
        assert report["ts_converged_reason"]>0
        assert report["ts_rejections"]==0 and report["snes_failures"]==0
        assert report["accepted_steps"]==4
        assert report["max_energy_formula_error"]<1e-8
        assert report["max_rank_local_full_state_discrepancy"]<1e-8
        assert trace["times"].shape==(5,)
        assert trace["states"].shape[0]==5 and trace["H"].shape==(5,)
        assert np.all(np.isfinite(trace["states"]))
        if baseline is None:
            baseline=trace
        else:
            assert np.array_equal(baseline["times"],trace["times"])
            assert baseline["states"].shape==trace["states"].shape
            whole=float(np.max(np.linalg.norm(
                baseline["states"]-trace["states"],axis=1
            )))
            Hdiff=float(np.max(abs(baseline["H"]-trace["H"])))
            if whole>=1e-8 or Hdiff>=1e-8:
                raise ValueError(f"Original example solver parity failed "
                                 f"{name}: whole_FE={whole}, H={Hdiff}")
            results[name]={"complete_FE_max_L2_against_default":whole,
                           "Hamiltonian_max_difference":Hdiff}
summary={
    "status":"PASS_ORIGINAL_LUMPED_EXAMPLES_MUMPS_SUPERLU_DIST_PARITY",
    "actual_native_cases":12,
    "examples": ["mass_spring_damper","nonlinear_pendulum"],
    "MPI_rank_counts":[1,2],"reference_backend":"mumps",
    "independent_comparisons_against_default":results,
    "published_git_snapshot_files_modified":False,
    "performance_benchmark_claimed":False,
    "full_default_simulation_horizon_qualified":False,
}
(root/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2))
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,allowZip64=True) as z:
    for p in sorted(root.rglob("*")):
        if p.is_file() and not any(
            x.startswith("work_") for x in p.relative_to(root).parts
        ):
            z.write(p,arcname=str(p.relative_to(root)))
print("SHAREABLE NATIVE ORIGINAL EXAMPLES ZIP:",archive)
PY
echo "COMPLETE original GetFEM/PETSc example LU pilot: $out"
