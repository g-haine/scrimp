#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Native diagnostic of published wave.py DAE initialization and complete
# 5-step CN trajectories, 1/2 MPI ranks. Does NOT alter the solver or
# rebuild/publish any Docker image. Output stays on the host, NEVER
# in the Docker build context. Dynamic data drift is measured, not hidden.
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then
  docker_cmd=(docker)
else
  docker_cmd=(sudo docker)
fi

image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
mesh="${SCRIMP_V130_DAE_MESH:-0.2}"
if [[ ! "$mesh" =~ ^0\.[0-9]+$ ]]; then
  echo "SCRIMP_V130_DAE_MESH must be 0.x (for example 0.2)" >&2
  exit 2
fi
for source in   sandbox/audit_v130_wave_dae_initial.py   sandbox/compare_v130_wave_dae_initial.py; do
  [[ -s "$source" ]] || { echo "Missing $source" >&2; exit 2; }
done

revision="$("${docker_cmd[@]}" image inspect --format   '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
if [[ -z "$revision" || "$revision" == "<no value>" ]]; then
  echo "Missing immutable OCI source revision label for $image" >&2
  exit 2
fi
head="$(git rev-parse HEAD)"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-dae-wave-h${mesh//./p}-${head:0:8}-XXXXXXXX")"
work="$(mktemp -d "$out/work-XXXXXXXX")"
echo "SCRIMP v1.3.0 physical DAE initial-state investigation"
echo "Image: $image (immutable source $revision)"
echo "Diagnostic script Git revision: $head"
echo "Published wave.py, h=$mesh, fixed CN dt=.002, final=.01, steps=5"
echo "MPI ranks: 1, 2; direct GetFEM->PETSc COO, exact MUMPS on_change"
echo "No production solver patch; no v1.2.0 rerun; no native compilation."
echo "Host evidence: $out"
printf '{"schema":"scrimp-v130-dae-run-v1","image":"%s","image_source_sha":"%s","audit_source_sha":"%s","mesh":%s,"steps":5,"dt":0.002,"t_final":0.01,"ranks":[1,2]}\n'   "$image" "$revision" "$head" "$mesh" > "$out/manifest.json"

for ranks in 1 2; do
  echo "=== RUN wave DAE audit: MPI=$ranks ==="
  if ! "${docker_cmd[@]}" run --rm     --user "$(id -u):$(id -g)"     --workdir /scrimp-work     --env HOME=/tmp --env MPLBACKEND=Agg     --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox     --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1     --env SCRIMP_IMAGE_REVISION="$revision"     --env SCRIMP_GETFEM_MATRIX_MODE=distributed     --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo     --env SCRIMP_MATRIX_MEMORY_MODE=compact     --env SCRIMP_TRAJECTORY_BACKEND=memory     --env SCRIMP_MUMPS_RESET_POLICY=on_change     --env SCRIMP_INIT_DYNAMIC_STATE_AUDIT=1     --env SCRIMP_INIT_FACTOR_SOLVER_TYPE=mumps     --volume "$work:/scrimp-work:rw"     --volume "$out:/artifacts:rw"     --volume "$root/sandbox/audit_v130_wave_dae_initial.py:/workspace/scrimp/sandbox/audit_v130_wave_dae_initial.py:ro"     "$image"     timeout -k 5s 240s       mpiexec --oversubscribe --bind-to none -n "$ranks"       python -u /workspace/scrimp/sandbox/audit_v130_wave_dae_initial.py         --mesh "$mesh" --dt 0.002 --steps 5 --final-time 0.01         --output-dir /artifacts     > "$out/wave_${ranks}r.log" 2>&1; then
    echo "FAILED: DAE audit MPI=$ranks, inspect $out/wave_${ranks}r.log" >&2
    tail -n 90 "$out/wave_${ranks}r.log" >&2
    exit 1
  fi
  grep -aE 'DIAGNOSTIC DAE|\\[SCRIMP initial dynamic audit\\]'     "$out/wave_${ranks}r.log" || true
done

echo "=== Compare actual complete time trajectories and DAE diagnostics ==="
if ! "${docker_cmd[@]}" run --rm   --user "$(id -u):$(id -g)"   --workdir /scrimp-work   --env HOME=/tmp --env PYTHONPATH=/workspace/scrimp/sandbox   --volume "$work:/scrimp-work:rw"   --volume "$out:/artifacts:rw"   --volume "$root/sandbox/compare_v130_wave_dae_initial.py:/workspace/scrimp/sandbox/compare_v130_wave_dae_initial.py:ro"   "$image"     python -u /workspace/scrimp/sandbox/compare_v130_wave_dae_initial.py       /artifacts   2>&1 | tee "$out/compare.log"; then
  echo "FAILED MPI parity or incomplete trajectory: $out" >&2
  exit 1
fi
echo "COMPLETE DIAGNOSTIC (not physical fidelity acceptance): $out"
