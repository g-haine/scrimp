#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# STRICT, SHORT A/B experiment: 3 CN steps, not a long simulation.
# One MPI Docker image, same physics, same output (memory backend),
# legacy versus compact PETSc matrix storage.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
sha="$(git rev-parse --short=8 HEAD)"
mkdir -p "$root/sandbox/outputs"
out="$(mktemp -d "$root/sandbox/outputs/compact-matrix-$sha-XXXXXX")"

echo "SCRIMP matrix memory compact-vs-legacy | commit=$sha"
echo "2 ranks, acoustic h=0.04 + nonlinear SWE h=0.35, 3 CN steps each"
echo "Output: $out"
bash -n sandbox/run_compact_matrix_parity.sh
sudo docker run --rm \
    --user "$(id -u):$(id -g)" \
    -e HOME=/tmp \
    -e MPLBACKEND=Agg \
    --cpuset-cpus="0,2" \
    -v "$root/scrimp:/workspace/scrimp/scrimp:ro" \
    -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
    -v "$root/examples:/workspace/scrimp/examples:ro" \
    -v "$root/tests:/workspace/scrimp/tests:ro" \
    -v "$out:/workspace/scrimp/outputs" \
    -v "$out:/artifacts" \
    scrimp-mpi:local \
    bash -c '
        set -euo pipefail
        echo "PREFLIGHT: pure-Python compact matrix Hessian algebra"
        python -m unittest -v tests.test_matrix_memory

        for case in acoustic swe; do
            for mode in legacy compact; do
                dest="/artifacts/$case/$mode"
                mkdir -p "$dest"
                echo "RUN: $case / $mode / 2 MPI ranks / 3 CN steps"
                if [[ "$case" == acoustic ]]; then
                    SCRIMP_MATRIX_MEMORY_MODE="$mode" \
                    SCRIMP_TRAJECTORY_BACKEND=memory \
                    SCRIMP_GETFEM_MATRIX_MODE=distributed \
                    SCRIMP_MUMPS_RESET_POLICY=never \
                    SCRIMP_BENCH_GIT_SHA="'"$sha"'" \
                    timeout -k 5s 180s mpiexec --bind-to none -n 2 \
                        python -u sandbox/verify_async_acoustic.py \
                            --mode memory \
                            --mesh 0.04 --dt 0.001 \
                            --final-time 0.003 \
                            --output-dir "$dest" \
                        > "$dest/solve.log" 2>&1 || {
                            tail -n 90 "$dest/solve.log" >&2
                            exit 1
                        }
                else
                    SCRIMP_MATRIX_MEMORY_MODE="$mode" \
                    SCRIMP_TRAJECTORY_BACKEND=memory \
                    SCRIMP_GETFEM_MATRIX_MODE=distributed \
                    SCRIMP_MUMPS_RESET_POLICY=on_change \
                    SCRIMP_BENCH_GIT_SHA="'"$sha"'" \
                    timeout -k 5s 240s mpiexec --bind-to none -n 2 \
                        python -u sandbox/verify_async_swe.py \
                            --mode memory \
                            --mesh 0.35 --dt 0.001 \
                            --final-time 0.003 \
                            --spinup-duration 0.03 \
                            --output-dir "$dest" \
                        > "$dest/solve.log" 2>&1 || {
                            tail -n 90 "$dest/solve.log" >&2
                            exit 1
                        }
                fi
                if [[ "$case" == acoustic ]]; then
                    suffix=physical
                else
                    suffix=swe
                fi
                test -s "$dest/memory_${suffix}.npz"
                test -s "$dest/memory_${suffix}.json"
                echo "  PASS: $case/$mode finished"
            done
        done

        echo "COMPARE: all coefficients, Hamiltonian, power and volume balances"
        python -u sandbox/compare_compact_matrix.py \
            --directory /artifacts
    ' || {
        echo "FAILED: inspect solve logs in $out/{acoustic,swe}/{legacy,compact}" >&2
        exit 1
    }

echo "SUCCESS: compact matrix-policy native parity gates completed."
echo "Results: $out/matrix_memory_comparison.json"
