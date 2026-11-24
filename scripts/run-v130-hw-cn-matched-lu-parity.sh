#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Four ORIGINAL Heat-Wave CN physical trajectories with strict boundary-
# compatible initial data: MUMPS vs experimental row-matched SuperLU_DIST.
# MPI1 only: original public ParaView export is not yet MPI2 qualified.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || {
 echo "Usage: $0" >&2;exit 2
}
docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker);fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
 --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
 echo "Unqualified native GetFEM/PETSc image: $revision" >&2;exit 2
}
python3 -m py_compile sandbox/compare_v130_hw_cn_matched_lu_parity.py
python3 tests/test_v130_hw_cn_matched_lu_parity.py
python3 tests/test_v130_hw_cn_matched_lu_parity.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-hw-cn-matched-parity-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-HW-original-CN-matched-parity-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "Original Heat-Wave full-state four-step CN, MPI1, strict algebraic DAE"
echo "Clamp and physical constrained-L2 p0; MUMPS vs row-matched SuperLU_DIST"
echo "Same GetFEM, COO, CN, Hamiltonians and real port powers"
echo "Scientific CN balance, solver parity AND non-vacuous interface power are separate gates"
echo "Published files immutable; no production row matching or TS policy change"
echo "EVIDENCE: $out"

for candidate in clamp constrained_l2; do
 for solver in mumps matched_superlu; do
  name="${candidate}_${solver}"
  mkdir "$out/$name" "$out/work_$name"
  if [[ "$solver" == "mumps" ]]; then
   target=/workspace/scrimp/sandbox/probe_v130_heat_wave_compatible_cn.py
   solver_args=(--candidate "$candidate" --backend mumps --out /artifacts)
  else
   target=/workspace/scrimp/sandbox/probe_v130_hw_matched_pc_cn.py
   solver_args=(--candidate "$candidate" --out /artifacts)
  fi
  echo "RUN original Heat-Wave candidate=$candidate solver=$solver"
  if ! "${docker_cmd[@]}" run --rm \
   --user "$(id -u):$(id -g)" --env HOME=/tmp --env MPLBACKEND=Agg \
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
   --volume "$root/examples/shallow_water.py:/workspace/scrimp/examples/shallow_water.py:ro" \
   --volume "$root/sandbox/probe_v130_main_2d_examples.py:/workspace/scrimp/sandbox/probe_v130_main_2d_examples.py:ro" \
   --volume "$root/sandbox/measure_v130_heat_wave_physics.py:/workspace/scrimp/sandbox/measure_v130_heat_wave_physics.py:ro" \
   --volume "$root/sandbox/probe_v130_heat_wave_boundary_compatibility.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_boundary_compatibility.py:ro" \
   --volume "$root/sandbox/probe_v130_heat_wave_compatible_cn.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_compatible_cn.py:ro" \
   --volume "$root/sandbox/probe_v130_hw_cn_identical_lu.py:/workspace/scrimp/sandbox/replay_base.py:ro" \
   --volume "$root/sandbox/probe_v130_hw_cn_exact_row_matching.py:/workspace/scrimp/sandbox/probe_v130_hw_cn_exact_row_matching.py:ro" \
   --volume "$root/sandbox/probe_v130_hw_matching_pc.py:/workspace/scrimp/sandbox/probe_v130_hw_matching_pc.py:ro" \
   --volume "$root/sandbox/probe_v130_hw_matched_pc_cn.py:/workspace/scrimp/sandbox/probe_v130_hw_matched_pc_cn.py:ro" \
   --volume "$out/$name:/artifacts:rw" \
   --volume "$out/work_$name:/scrimp-work:rw" --workdir /scrimp-work \
   "$image" timeout -k 5s "${SCRIMP_V130_HW_MATCHED_PARITY_TIMEOUT:-300s}" \
   mpiexec --oversubscribe --bind-to none -n 1 \
   python -u "$target" "${solver_args[@]}" \
   > "$out/$name.host.log" 2>&1; then
   echo "FAIL actual CN run $name; original outputs retained: $out" >&2
   tail -n 160 "$out/$name.host.log" >&2
   exit 1
  fi
 done
done
python3 sandbox/compare_v130_hw_cn_matched_lu_parity.py --root "$out"
echo "COMPLETE four-way actual Heat-Wave CN backend parity: $out"
