#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp
#
# Very short native-MPI transport smoke: no PDE or long time simulation.
# The writer is launched by the container SHELL before mpiexec starts,
# so the MPI workers never spawn/fork an auxiliary process.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 1
root="$PWD"
sha="$(git rev-parse --short=8 HEAD)"
out="$root/sandbox/outputs/async-binary-smoke-$sha"
mkdir -p "$out"

bash -n sandbox/run_async_binary_smoke.sh
echo "SCRIMP asynchronous binary I/O smoke | commit=$sha | MPI=2 | four exact snapshots"
echo "Output: $out"

sudo docker run --rm \
    --user "$(id -u):$(id -g)" \
    -e HOME=/tmp \
    --cpuset-cpus="0,2" \
    -v "$root/sandbox:/workspace/scrimp/sandbox:ro" \
    -v "$root/tests:/workspace/scrimp/tests:ro" \
    -v "$out:/artifacts" \
    scrimp-mpi:local \
    bash -c '
        set -euo pipefail
        echo "PREFLIGHT: subprocess binary format, finite queue, parity tests"
        python -m unittest -v tests.test_async_binary_trajectory

        echo "WRITER: starting external process before mpiexec"
        python -u sandbox/async_binary_trajectory.py \
            --writer --socket /artifacts/trajectory_writer.sock \
            --output-dir /artifacts --ranks 2 --timeout 30 \
            > /artifacts/writer.log 2>&1 &
        writer_pid=$!
        trap '\''kill "$writer_pid" 2>/dev/null || true'\'' EXIT

        echo "MPI: stream only each owned PETSc Vec slice"
        mpiexec --bind-to none -n 2 \
            python -u sandbox/probe_async_binary_mpi.py \
            --socket /artifacts/trajectory_writer.sock --n-global 17

        wait "$writer_pid"
        trap - EXIT
        echo "READBACK: verify every coefficient, timestep and partition"
        python sandbox/check_async_binary_mpi.py \
            --directory /artifacts --ranks 2 --n-global 17
    '

echo "PASS: isolated asynchronous partitioned binary transport."
echo "Files: $out"
