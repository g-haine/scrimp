#!/usr/bin/env bash
# SCRIMP 1.3: original Heat-Wave CN, FULL 15-second MPI1 vs MPI2 parity.
# This gate is DISTINCT from the already qualified four-step MPI parity
# and MPI1 3000-step physical power balance (source 97c7606b).
# Preserve EVERY original GetFEM FE state for complete 3001x7689 parity.
# No recomputation of already accepted VTK/PNG/GIF; no parallel PV claim.
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
    echo "Unqualified native GetFEM/PETSc image revision: $revision" >&2
    exit 2
}

python3 -m py_compile \
    examples/heat_wave.py \
    sandbox/measure_v130_hw_dual_facet_physics.py \
    sandbox/probe_v130_heat_wave_mpi_extended.py \
    scripts/verify-v130-heat-wave-full-mpi2.py \
    scripts/reassess-v130-heat-wave-full-mpi2-norm-parity.py \
    scripts/verify-v130-heat-wave-full-mpi2-normwise.py \
    tests/test_v130_heat_wave_full_mpi2_source.py \
    tests/test_v130_heat_wave_full_mpi_norm_release_source.py
python3 tests/test_v130_heat_wave_full_mpi2_source.py
python3 tests/test_v130_heat_wave_full_mpi_norm_release_source.py
python3 tests/test_v130_heat_wave_mpi_long_source.py
python3 scripts/check-v130-published-immutable.py
bash -n "$0"

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-heat-wave-full-15s-MPI1-MPI2-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
python3 - "$out" "$image" "$revision" <<'PY'
import json,subprocess,sys
from pathlib import Path
root=Path(sys.argv[1])
(root/"manifest.json").write_text(json.dumps({
    "schema":"scrimp-v130-heat-wave-full-MPI1-MPI2-native-15s-v1",
    "source":subprocess.check_output(
        ["git","rev-parse","HEAD"],text=True
    ).strip(),
    "image":sys.argv[2],
    "image_revision":sys.argv[3],
    "canonical_example":"examples/heat_wave.py",
    "original_model_unchanged":True,
    "original_conforming_two_sided_interface":True,
    "physics_only_no_PV_writer":True,
    "profile":"gaussian",
    "dt":.005,
    "steps":3000,
    "tf":15.,
    "FEM_DOFs":7689,
    "full_FE_states_saved_per_rank":True,
    "MPI":[2,1],
    "previously_qualified_short_MPI1_MPI2":True,
    "previously_qualified_MPI1_long_power":True,
    "current_MPI2_full_qualified":False,
    "parallel_VTK_qualified":False,
    "other_solver_qualified":False,
    "not_a_modal_eigenmode_proof":True,
},indent=2)+"\n")
PY
echo "Heat-Wave original full fifteen-second MUMPS MPI1 / MPI2 physical parity"
echo "Real GetFEM 7689 DOFs, 3000 fixed CN steps dt=.005, 3001 genuine FE states"
echo "FIRST run MPI2 (new scientific question); THEN MPI1 (independent baseline)"
echo "All original GetFEM power observables, strict DAE and full coefficient VECTOR normwise parity"
echo "Archiving all 3x3001x7689 original FE coefficients; evidence may be large"
echo "No synthetic FE states, no parallel ParaView output, no gallery-mode claim"
echo "EVIDENCE: $out"

run_case() {
    local label="$1" ranks="$2"
    mkdir -p "$out/$label" "$out/work_$label"
    echo "RUN canonical Heat-Wave CN full: $label, ranks=$ranks"
    local rc=0
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
        --volume "$root/examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro" \
        --volume "$root/sandbox/measure_v130_hw_dual_facet_physics.py:/workspace/scrimp/sandbox/measure_v130_hw_dual_facet_physics.py:ro" \
        --volume "$root/sandbox/probe_v130_heat_wave_mpi_extended.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_mpi_extended.py:ro" \
        --volume "$out/$label:/artifacts:rw" \
        --volume "$out/work_$label:/scrimp-work:rw" \
        --workdir /scrimp-work "$image" \
        timeout -k 10s "${SCRIMP_V130_HW_FULL_MPI2_TIMEOUT:-7200s}" \
        mpiexec --oversubscribe --bind-to none -n "$ranks" \
        python -u /workspace/scrimp/sandbox/probe_v130_heat_wave_mpi_extended.py \
            --out /artifacts --dt 0.005 --steps 3000 \
            --profile gaussian --store-full-states \
        >"$out/$label.host.log" 2>&1 || rc=$?
    printf '%s\n' "$rc" > "$out/$label.exitcode"
    if [[ "$rc" != 0 ]]; then
        echo "UNQUALIFIED real $label exit=$rc; native files retained" >&2
        tail -n 85 "$out/$label.host.log" >&2
    else
        echo "NATIVE COMPLETE original full $label"
    fi
    return "$rc"
}

overall=0
if run_case mpi2_long 2; then
    if ! run_case mpi1_long 1; then overall=1; fi
else
    overall=1
    echo "MPI2 long failed; already-qualified MPI1 long not rerun" >&2
    echo "SKIPPED_AFTER_MPI2_FAILURE" > "$out/mpi1_long.exitcode"
fi

compare_rc=0
# The historical strict per-DOF checker is preserved as a separate
# scientific diagnostic, NOT the current release criterion.
# The current reusable gate independently checks genuine native
# physical power traces and all 3001 original FE coefficient vectors
# using the previously accepted normwise parity metric.
python3 scripts/verify-v130-heat-wave-full-mpi2-normwise.py "$out" \
    --out "$out/summary.json" || compare_rc=$?
if [[ "$compare_rc" != 0 ]]; then overall=1; fi

# Preserve the complete original native per-rank FE vectors despite
# archive size. In particular, archive MPI2 rank0 AND rank1.
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1])
summary=json.loads((root/"summary.json").read_text())
archive=Path(str(root)+".zip")
with zipfile.ZipFile(
    archive,"x",compression=zipfile.ZIP_DEFLATED,
    compresslevel=2,allowZip64=True,
) as z:
    for f in sorted(root.rglob("*")):
        if f.is_file():
            z.write(f,arcname=str(f.relative_to(root)))
with zipfile.ZipFile(archive) as z:
    faulty=z.testzip()
    if faulty is not None:
        raise RuntimeError("Heat-Wave original MPI2 complete ZIP CRC failed: "+faulty)
print(json.dumps({
    "status":summary["status"],
    "errors":summary["errors"],
    "MPI2_full_horizon_qualified":summary["MPI2_full_horizon_qualified"],
    "max_FEM_coefficient_difference":summary.get("max_FEM_coefficient_difference"),
    "evidence":str(root),
    "archive":str(archive),
    "ZIP_CRC_pass":True,
},indent=2),flush=True)
PY
if [[ "$overall" == 0 ]]; then
    echo "PASS_HEAT_WAVE_FULL_15S_MPI1_MPI2_CN_NORM_PARITY"
else
    echo "FAIL_OR_UNQUALIFIED_FULL_HEAT_WAVE_MPI2 — archive preserved" >&2
fi
exit "$overall"
