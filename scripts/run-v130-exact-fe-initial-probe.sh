#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# SCIENTIFIC POSITIVE CONTROL ONLY: published wave p0 replaced in a
# runtime-only interception by p0=x*(3-x), e_p=x, q0=0, U_L(0)=0.
# Test whether exactly FE-representable compatible initial data admit
# an algebraic-only correction with every prescribed q/p DOF frozen.
# Never touch production code or rebuild the qualified rc3 image.
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
for src in \
  sandbox/probe_v130_wave_fixed_dynamics.py \
  sandbox/solve_v130_dae_algebraic_block.py \
  sandbox/audit_v130_wave_dae_initial.py \
  sandbox/compare_v130_wave_fixed_dynamics.py; do
  [[ -s "$src" ]] || { echo "Missing $src" >&2; exit 2; }
done
revision="$("${docker_cmd[@]}" image inspect --format \
  '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
expected="1d6d0c90115213298047b28237a7f38087563a35"
if [[ "$revision" != "$expected" ]]; then
  echo "This probe requires the exact original rc3 source: $expected" >&2
  echo "Found $revision on image $image" >&2
  exit 2
fi
head="$(git rev-parse HEAD)"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-dae-fe-positive-h${mesh//./p}-${head:0:8}-XXXXXXXX")"
work="$(mktemp -d "$out/work-XXXXXXXX")"
echo "SCRIMP v1.3.0 DAE FE-compatible positive control"
echo "Image $image @ $revision; Git diagnostic $head"
echo "Only p0 replaced: p0=x*(3-x), e_p=x, left U_L(0)=0"
echo "All original bricks, boundary controls, FE spaces, CN and MUMPS unchanged"
echo "MPI 1/2, h=$mesh, dt=0.002, 5 CN steps"
echo "Output: $out"
printf '{"schema":"scrimp-v130-dae-fe-positive-v1","source":"%s","diagnostic":"%s","mesh":"%s","scenario":"polynomial_compatible"}\n' \
  "$revision" "$head" "$mesh" > "$out/manifest.json"
common=( \
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
)
echo "PREFLIGHT: build exact-FE initial momentum in GetFEM"
if ! "${docker_cmd[@]}" run --rm "${common[@]}" "$image" \
  python -u /workspace/scrimp/sandbox/probe_v130_wave_fixed_dynamics.py \
    --mesh "$mesh" --scenario polynomial_compatible \
    --output-dir /artifacts --preflight-control-only \
  > "$out/preflight.log" 2>&1; then
  echo "FAIL: GetFEM positive-control construction: $out/preflight.log" >&2
  tail -n 90 "$out/preflight.log" >&2
  exit 1
fi
grep -aF "PASS: GetFEM constructed scenario=polynomial_compatible" "$out/preflight.log"
for ranks in 1 2; do
  echo "RUN native FE-positive DAE MPI=$ranks"
  if ! "${docker_cmd[@]}" run --rm "${common[@]}" "$image" \
      timeout -k 5s 240s \
      mpiexec --oversubscribe --bind-to none -n "$ranks" \
      python -u /workspace/scrimp/sandbox/probe_v130_wave_fixed_dynamics.py \
        --mesh "$mesh" --scenario polynomial_compatible \
        --output-dir /artifacts \
      > "$out/fe_positive_${ranks}r.log" 2>&1; then
    echo "FAIL: FE-positive DAE MPI=$ranks; evidence $out" >&2
    tail -n 90 "$out/fe_positive_${ranks}r.log" >&2
    exit 1
  fi
  grep -aF "DAE ALGEBRAIC-ONLY" "$out/fe_positive_${ranks}r.log" || true
done
echo "COMPARE: fixed q/p exact DAE feasibility and pseudo dynamic drift"
if ! "${docker_cmd[@]}" run --rm \
    --user "$(id -u):$(id -g)" \
    --workdir /scrimp-work \
    --env HOME=/tmp --env PYTHONPATH=/workspace/scrimp/sandbox \
    --volume "$work:/scrimp-work:rw" \
    --volume "$out:/artifacts:rw" \
    --volume "$root/sandbox/compare_v130_wave_fixed_dynamics.py:/workspace/scrimp/sandbox/compare_v130_wave_fixed_dynamics.py:ro" \
    "$image" \
      python -u /workspace/scrimp/sandbox/compare_v130_wave_fixed_dynamics.py \
        /artifacts --scenario polynomial_compatible \
    2>&1 | tee "$out/compare_fe_positive.log"; then
  echo "FAIL: incompatible MPI feasibility conclusions; evidence $out" >&2
  exit 1
fi
echo "COMPLETE: FE-compatible native DAE positive control; output $out"
