#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Tiny diagnostic-only full CN and Hamiltonian/port-power A/B.
# ORIGINAL qualified rc3 image, no production code/image changes.
# Run SAME nonzero compatible initial data and driven boundary with
# pseudo-TS state or algebraic-only projected main CN start at MPI 1/2.
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
  echo "Expected SCRIMP_V130_DAE_MESH=0.x" >&2
  exit 2
fi
for file in \
  sandbox/audit_v130_wave_initializers.py \
  sandbox/compare_v130_wave_initializers.py \
  sandbox/probe_v130_wave_fixed_dynamics.py \
  sandbox/solve_v130_dae_algebraic_block.py \
  sandbox/audit_v130_wave_dae_initial.py \
  sandbox/check_acoustic_energy_mpi.py; do
  [[ -s "$file" ]] || { echo "Missing $file" >&2; exit 2; }
done
revision="$("${docker_cmd[@]}" image inspect --format \
  '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
expected="1d6d0c90115213298047b28237a7f38087563a35"
if [[ "$revision" != "$expected" ]]; then
  echo "Refusing unqualified image revision $revision (expected $expected)" >&2
  exit 2
fi
head="$(git rev-parse HEAD)"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-wave-init-power-h${mesh//./p}-${head:0:8}-XXXXXXXX")"
work="$(mktemp -d "$out/work-XXXXXXXX")"
echo "SCRIMP v1.3.0 | initial DAE projection versus pseudo-TS | CN + physical port work"
echo "Native runtime image: $image @ $revision"
echo "Diagnostic Git revision: $head"
echo "Same h=$mesh, p0=x*(3-x), q0=0, U_L=0.08*sin(20*t)"
echo "dt=.002, t_f=.02, 10 fixed CN steps, MPI 1/2"
echo "4 runs: pseudo/algebraic x 1/2 ranks, full states, Hamiltonian and signed port power"
echo "No production changes, no image rebuild, no expensive scaling rerun"
echo "Host evidence: $out"
printf '{"schema":"scrimp-v130-wave-initializer-physics-run-v1","image_revision":"%s","diagnostic_revision":"%s","mesh":"%s","initial_p0":"x*(3-x)","drive_left":"0.08*sin(20*t)","steps":10,"dt":0.002,"mpi":[1,2]}\n' \
  "$revision" "$head" "$mesh" > "$out/manifest.json"
common=(\
  --user "$(id -u):$(id -g)" \
  --workdir /scrimp-work \
  --env HOME=/tmp --env MPLBACKEND=Agg \
  --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
  --env SCRIMP_IMAGE_REVISION="$revision" \
  --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
  --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
  --env SCRIMP_MATRIX_MEMORY_MODE=compact \
  --env SCRIMP_TRAJECTORY_BACKEND=memory \
  --env SCRIMP_MUMPS_RESET_POLICY=on_change \
  --env SCRIMP_INIT_FACTOR_SOLVER_TYPE=mumps \
  --volume "$work:/scrimp-work:rw" \
  --volume "$out:/artifacts:rw" \
  --volume "$root/sandbox/audit_v130_wave_initializers.py:/workspace/scrimp/sandbox/audit_v130_wave_initializers.py:ro" \
  --volume "$root/sandbox/probe_v130_wave_fixed_dynamics.py:/workspace/scrimp/sandbox/probe_v130_wave_fixed_dynamics.py:ro" \
  --volume "$root/sandbox/audit_v130_wave_dae_initial.py:/workspace/scrimp/sandbox/audit_v130_wave_dae_initial.py:ro" \
  --volume "$root/sandbox/solve_v130_dae_algebraic_block.py:/workspace/scrimp/sandbox/solve_v130_dae_algebraic_block.py:ro" \
  --volume "$root/sandbox/check_acoustic_energy_mpi.py:/workspace/scrimp/sandbox/check_acoustic_energy_mpi.py:ro" \
)
for mode in pseudo algebraic; do
  for ranks in 1 2; do
    echo "RUN mode=$mode ranks=$ranks"
    if ! "${docker_cmd[@]}" run --rm "${common[@]}" "$image" \
      timeout -k 5s 240s \
      mpiexec --oversubscribe --bind-to none -n "$ranks" \
      python -u /workspace/scrimp/sandbox/audit_v130_wave_initializers.py \
        --mode "$mode" --mesh "$mesh" \
        --dt 0.002 --steps 10 \
        --drive-left "0.08*sin(20*t)" \
        --output-dir /artifacts \
      > "$out/${mode}_${ranks}r.log" 2>&1; then
      echo "FAIL: native initializer $mode/$ranks. Inspect $out/${mode}_${ranks}r.log" >&2
      tail -n 90 "$out/${mode}_${ranks}r.log" >&2
      exit 1
    fi
    grep -aF "SCIENCE $mode MPI=$ranks" "$out/${mode}_${ranks}r.log" || true
  done
done
echo "COMPARE full trajectories and physical midpoint port powers at MPI 1/2"
if ! "${docker_cmd[@]}" run --rm \
    --user "$(id -u):$(id -g)" \
    --workdir /scrimp-work \
    --env HOME=/tmp --env PYTHONPATH=/workspace/scrimp/sandbox \
    --volume "$work:/scrimp-work:rw" \
    --volume "$out:/artifacts:rw" \
    --volume "$root/sandbox/compare_v130_wave_initializers.py:/workspace/scrimp/sandbox/compare_v130_wave_initializers.py:ro" \
    "$image" \
      python -u /workspace/scrimp/sandbox/compare_v130_wave_initializers.py \
        /artifacts \
    2>&1 | tee "$out/compare.log"; then
  echo "FAIL: scientific initializer comparison; inspect $out/compare.log" >&2
  exit 1
fi
echo "COMPLETE: physical CN initialization A/B; host evidence $out"
