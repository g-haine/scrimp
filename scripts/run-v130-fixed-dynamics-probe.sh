#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Numerical feasibility of correcting ONLY algebraic wave variables at t0,
# with prescribed dynamic q,p coefficients frozen; original rc3 solver and
# accepted pseudo-TS initial state are NEVER changed.
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
for source in \
    sandbox/probe_v130_wave_fixed_dynamics.py \
    sandbox/solve_v130_dae_algebraic_block.py \
    sandbox/audit_v130_wave_dae_initial.py \
    sandbox/compare_v130_wave_fixed_dynamics.py \
    sandbox/compare_v130_wave_control_compatibility.py; do
  [[ -s "$source" ]] || { echo "Missing $source" >&2; exit 2; }
done
revision="$("${docker_cmd[@]}" image inspect --format \
    '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
expected="1d6d0c90115213298047b28237a7f38087563a35"
if [[ "$revision" != "$expected" ]]; then
  echo "Only exact qualified rc3 image is allowed: $revision" >&2
  exit 2
fi
head="$(git rev-parse HEAD)"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-algebraic-only-h${mesh//./p}-${head:0:8}-XXXXXXXX")"
work="$(mktemp -d "$out/work-XXXXXXXX")"
echo "SCRIMP v1.3.0 initial DAE algebraic-only feasibility probe"
echo "Docker image: $image @ $revision; diagnostic Git source $head"
echo "Fixed differential q(0), p(0); linear wave IFunction algebraic rows"
echo "Same published wave, h=$mesh, CN dt=0.002 t_f=0.01, MPI 1/2"
echo "Numeric feasibility + actual PETSc IFunction re-evaluation"
echo "No production change, no solver rewrite, no Docker rebuild."
echo "Host output: $out"
printf '{"schema":"scrimp-v130-algebraic-only-run-v1","image_revision":"%s","diagnostic_revision":"%s","mesh":"%s"}\n' \
  "$revision" "$head" "$mesh" > "$out/manifest.json"
# Fail early on GetFEM grammar errors in the diagnostic control.
# The earlier A/B run reached published 1/2 ranks then stopped while
# parsing the matched source-brick expression; avoid repeating that.
echo "PREFLIGHT: GetFEM parses matched-left control source (MPI=1)"
if ! "${docker_cmd[@]}" run --rm \
  --user "$(id -u):$(id -g)" \
  --workdir /scrimp-work \
  --env HOME=/tmp --env MPLBACKEND=Agg \
  --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
  --volume "$work:/scrimp-work:rw" \
  --volume "$out:/artifacts:rw" \
  --volume "$root/sandbox/probe_v130_wave_fixed_dynamics.py:/workspace/scrimp/sandbox/probe_v130_wave_fixed_dynamics.py:ro" \
  --volume "$root/sandbox/solve_v130_dae_algebraic_block.py:/workspace/scrimp/sandbox/solve_v130_dae_algebraic_block.py:ro" \
  --volume "$root/sandbox/audit_v130_wave_dae_initial.py:/workspace/scrimp/sandbox/audit_v130_wave_dae_initial.py:ro" \
  "$image" \
    python -u /workspace/scrimp/sandbox/probe_v130_wave_fixed_dynamics.py \
      --mesh "$mesh" --scenario matched_left_control \
      --output-dir /artifacts --preflight-control-only \
  > "$out/control_preflight.log" 2>&1; then
  echo "ERROR: matched-left GetFEM control source did not parse." >&2
  tail -n 90 "$out/control_preflight.log" >&2
  echo "Host preflight log: $out/control_preflight.log" >&2
  exit 1
fi
grep -aF "PASS: GetFEM parsed the matched-left" \
  "$out/control_preflight.log"

for scenario in published matched_left_control; do
  for ranks in 1 2; do
  label="${ranks}r_${scenario}"
  echo "RUN projection feasibility MPI=$ranks scenario=$scenario"
  if ! "${docker_cmd[@]}" run --rm \
    --user "$(id -u):$(id -g)" \
    --workdir /scrimp-work \
    --env HOME=/tmp --env MPLBACKEND=Agg \
    --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
    --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
    --env SCRIMP_IMAGE_REVISION="$revision" \
    --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
    --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
    --env SCRIMP_MATRIX_MEMORY_MODE=compact \
    --env SCRIMP_MUMPS_RESET_POLICY=on_change \
    --env SCRIMP_INIT_FACTOR_SOLVER_TYPE=mumps \
    --volume "$work:/scrimp-work:rw" \
    --volume "$out:/artifacts:rw" \
    --volume "$root/sandbox/probe_v130_wave_fixed_dynamics.py:/workspace/scrimp/sandbox/probe_v130_wave_fixed_dynamics.py:ro" \
    --volume "$root/sandbox/solve_v130_dae_algebraic_block.py:/workspace/scrimp/sandbox/solve_v130_dae_algebraic_block.py:ro" \
    --volume "$root/sandbox/audit_v130_wave_dae_initial.py:/workspace/scrimp/sandbox/audit_v130_wave_dae_initial.py:ro" \
    "$image" \
    timeout -k 5s 240s \
      mpiexec --oversubscribe --bind-to none -n "$ranks" \
      python -u /workspace/scrimp/sandbox/probe_v130_wave_fixed_dynamics.py \
        --mesh "$mesh" --output-dir /artifacts \
        --scenario "$scenario" \
    > "$out/probe_${label}.log" 2>&1; then
    echo "ERROR: projection feasibility run failed MPI=$ranks scenario=$scenario: $out" >&2
    tail -n 90 "$out/probe_${label}.log" >&2
    exit 1
  fi
  grep -aF "DAE ALGEBRAIC-ONLY" "$out/probe_${label}.log" || true
  done
  echo "COMPARE 1/2 MPI feasibility for scenario=$scenario"
if ! "${docker_cmd[@]}" run --rm \
  --user "$(id -u):$(id -g)" \
  --workdir /scrimp-work \
  --env HOME=/tmp --env PYTHONPATH=/workspace/scrimp/sandbox \
  --volume "$work:/scrimp-work:rw" \
  --volume "$out:/artifacts:rw" \
  --volume "$root/sandbox/compare_v130_wave_fixed_dynamics.py:/workspace/scrimp/sandbox/compare_v130_wave_fixed_dynamics.py:ro" \
  "$image" \
    python -u /workspace/scrimp/sandbox/compare_v130_wave_fixed_dynamics.py \
      /artifacts --scenario "$scenario" \
  2>&1 | tee "$out/compare_${scenario}.log"; then
  echo "ERROR: algebraic feasibility comparison failed: $out" >&2
  exit 1
  fi
done

echo "COMPARE original published control vs trace-matched left boundary"
if ! "${docker_cmd[@]}" run --rm \
  --user "$(id -u):$(id -g)" \
  --workdir /scrimp-work \
  --env HOME=/tmp --env PYTHONPATH=/workspace/scrimp/sandbox \
  --volume "$work:/scrimp-work:rw" \
  --volume "$out:/artifacts:rw" \
  --volume "$root/sandbox/compare_v130_wave_control_compatibility.py:/workspace/scrimp/sandbox/compare_v130_wave_control_compatibility.py:ro" \
  "$image" \
    python -u /workspace/scrimp/sandbox/compare_v130_wave_control_compatibility.py \
      /artifacts \
  2>&1 | tee "$out/compare_control_compatibility.log"; then
  echo "ERROR: left-control compatibility contrast failed: $out" >&2
  exit 1
fi
echo "COMPLETE: interpretation-dependent feasibility diagnostic, logs at $out"
