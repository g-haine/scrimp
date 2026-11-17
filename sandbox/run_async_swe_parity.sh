#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# SHORT physical test, not a long simulation:
# 2 MPI ranks, h=.35, 3 fixed CN steps, 2x same nonlinear SWE spin-up model.
# Independent writer started BEFORE mpiexec, never from within a rank.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
sha="$(git rev-parse --short=8 HEAD)"
mkdir -p "$root/sandbox/outputs"
out="$(mktemp -d "$root/sandbox/outputs/async-swe-parity-$sha-XXXXXX")"

echo "SCRIMP short nonlinear Saint-Venant async trajectory parity | commit=$sha"
echo "2 MPI ranks, h=0.35, dt=0.001, t_f=0.003 (3 CN steps, 0.03s spin-up)"
echo "Output: $out"
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

        echo "PREFLIGHT: independent binary writer protocol tests"
        python -m unittest -q tests.test_async_binary_trajectory

        echo "PREFLIGHT: load complete Saint-Venant driver before mpiexec"
        if ! python -u sandbox/verify_async_swe.py --help \
            > /artifacts/driver_import_preflight.log 2>&1; then
            echo "FAILED: Saint-Venant driver import; no MPI started"
            cat /artifacts/driver_import_preflight.log
            exit 1
        fi
        echo "  PASS: Saint-Venant model and diagnostics imports"

        echo "REFERENCE: PETSc + in-RAM snapshots"
        env SCRIMP_GETFEM_MATRIX_MODE=distributed \
            SCRIMP_MUMPS_RESET_POLICY=on_change \
            SCRIMP_TRAJECTORY_BACKEND=memory \
            SCRIMP_BENCH_GIT_SHA="'"$sha"'" \
        timeout -k 5s 300s mpiexec --bind-to none -n 2 \
            python -u sandbox/verify_async_swe.py \
            --mode memory --output-dir /artifacts \
            --mesh 0.35 --dt 0.001 --final-time 0.003 \
            --spinup-duration 0.03 \
            > /artifacts/memory.log 2>&1 || {
                tail -n 70 /artifacts/memory.log
                exit 1
            }

        echo "STREAMING: start separate I/O process before mpiexec"
        mkdir -p /artifacts/binary
        python -u sandbox/async_binary_trajectory.py \
            --writer --socket /artifacts/binary/trajectory_writer.sock \
            --output-dir /artifacts/binary --ranks 2 --timeout 300 \
            > /artifacts/binary/writer.log 2>&1 &
        writer_pid=$!
        trap '\''kill "$writer_pid" 2>/dev/null || true'\'' EXIT

        env SCRIMP_GETFEM_MATRIX_MODE=distributed \
            SCRIMP_MUMPS_RESET_POLICY=on_change \
            SCRIMP_TRAJECTORY_BACKEND=async_binary \
            SCRIMP_TRAJECTORY_DIR=/artifacts/binary \
            SCRIMP_TRAJECTORY_SOCKET=/artifacts/binary/trajectory_writer.sock \
            SCRIMP_TRAJECTORY_MAX_QUEUED=8 \
            SCRIMP_BENCH_GIT_SHA="'"$sha"'" \
        timeout -k 5s 300s mpiexec --bind-to none -n 2 \
            python -u sandbox/verify_async_swe.py \
            --mode async_binary --output-dir /artifacts \
            --mesh 0.35 --dt 0.001 --final-time 0.003 \
            --spinup-duration 0.03 \
            > /artifacts/async_binary.log 2>&1 || {
                tail -n 70 /artifacts/async_binary.log
                cat /artifacts/binary/writer.log
                exit 1
            }

        wait "$writer_pid" || {
            echo "External writer failed to commit"
            cat /artifacts/binary/writer.log
            exit 1
        }
        trap - EXIT
        echo "NUMERICAL: compare ALL states, Hamiltonian and computed port powers"
        python sandbox/compare_async_swe.py --directory /artifacts
    ' || {
        echo "FAILED: consult $out/{memory,async_binary}.log and binary/writer.log" >&2
        exit 1
    }
echo "PASS: short nonlinear Saint-Venant physical parity with independent I/O process"
echo "Output: $out"
