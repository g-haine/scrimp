#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Diagnostic-only, native distributed PETSc LSQR DAE projection.
# Original qualified rc3 runtime, no Docker rebuild and NO pseudo-TS.
# Positive polynomial: full CN + nonzero boundary port energy check.
# Negative published Gaussian: expect an explicit t0 rejection.
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
  echo "SCRIMP_V130_DAE_MESH must be 0.x" >&2
  exit 2
fi
for source in \
    sandbox/probe_v130_petsc_mpi_initial.py \
    sandbox/compare_v130_petsc_mpi_initial.py \
    sandbox/probe_v130_wave_fixed_dynamics.py \
    sandbox/audit_v130_wave_dae_initial.py \
    sandbox/audit_v130_wave_initializers.py \
    sandbox/solve_v130_dae_algebraic_block.py \
    sandbox/check_acoustic_energy_mpi.py; do
  [[ -s "$source" ]] || { echo "Missing source $source" >&2; exit 2; }
done
revision="$("${docker_cmd[@]}" image inspect --format \
    '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
expected="1d6d0c90115213298047b28237a7f38087563a35"
if [[ "$revision" != "$expected" ]]; then
  echo "Refusing unqualified Docker image; source=$revision" >&2
  exit 2
fi
head="$(git rev-parse HEAD)"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-petsc-dae-mpi-h${mesh//./p}-${head:0:8}-XXXXXXXX")"
work="$(mktemp -d "$out/work-XXXXXXXX")"
echo "SCRIMP v1.3.0 fully PETSc-distributed DAE initial feasibility"
echo "Qualified image $image source $revision; diagnostic Git $head"
echo "MPI=1/2; positive polynomial and incompatible Gaussian; no pseudo-TS"
echo "10 fixed CN steps and signed port work for positive data"
echo "Original source and Docker image never altered"
echo "Host evidence: $out"
printf '{"schema":"scrimp-v130-petsc-distributed-dae-v1","image_revision":"%s","diagnostic_revision":"%s","mesh":"%s","mpi":[1,2],"scenarios":["polynomial_compatible","published"]}\n' \
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
    --volume "$root/sandbox/probe_v130_petsc_mpi_initial.py:/workspace/scrimp/sandbox/probe_v130_petsc_mpi_initial.py:ro" \
    --volume "$root/sandbox/probe_v130_wave_fixed_dynamics.py:/workspace/scrimp/sandbox/probe_v130_wave_fixed_dynamics.py:ro" \
    --volume "$root/sandbox/audit_v130_wave_dae_initial.py:/workspace/scrimp/sandbox/audit_v130_wave_dae_initial.py:ro" \
    --volume "$root/sandbox/audit_v130_wave_initializers.py:/workspace/scrimp/sandbox/audit_v130_wave_initializers.py:ro" \
    --volume "$root/sandbox/solve_v130_dae_algebraic_block.py:/workspace/scrimp/sandbox/solve_v130_dae_algebraic_block.py:ro" \
    --volume "$root/sandbox/check_acoustic_energy_mpi.py:/workspace/scrimp/sandbox/check_acoustic_energy_mpi.py:ro" \
)
for scenario in polynomial_compatible published; do
  for ranks in 1 2; do
    echo "RUN PETSc fully distributed: $scenario MPI=$ranks"
    if ! "${docker_cmd[@]}" run --rm "${common[@]}" "$image" \
        timeout -k 5s 240s \
        mpiexec --oversubscribe --bind-to none -n "$ranks" \
        python -u /workspace/scrimp/sandbox/probe_v130_petsc_mpi_initial.py \
          --scenario "$scenario" --mesh "$mesh" \
          --dt 0.002 --steps 10 --output-dir /artifacts \
        > "$out/${scenario}_${ranks}r.log" 2>&1; then
      echo "ERROR: PETSc distributed DAE case failed: $scenario MPI=$ranks" >&2
      tail -n 90 "$out/${scenario}_${ranks}r.log" >&2
      echo "Host evidence: $out" >&2
      exit 1
    fi
    grep -aF "DISTRIBUTED PETSC DAE MPI=$ranks" "$out/${scenario}_${ranks}r.log" || true
  done
done
echo "COMPARE strict MPI/full-CN/port-power positive case and Gaussian rejection"
if ! "${docker_cmd[@]}" run --rm \
    --user "$(id -u):$(id -g)" \
    --workdir /scrimp-work \
    --env HOME=/tmp --env PYTHONPATH=/workspace/scrimp/sandbox \
    --volume "$work:/scrimp-work:rw" \
    --volume "$out:/artifacts:rw" \
    --volume "$root/sandbox/compare_v130_petsc_mpi_initial.py:/workspace/scrimp/sandbox/compare_v130_petsc_mpi_initial.py:ro" \
    "$image" \
      python -u /workspace/scrimp/sandbox/compare_v130_petsc_mpi_initial.py \
        /artifacts \
    2>&1 | tee "$out/compare.log"; then
  echo "ERROR: strict PETSc distributed DAE acceptance failed; see $out" >&2
  exit 1
fi
echo "COMPLETE: PETSc MPI-only DAE projection scientific probe; output $out"
