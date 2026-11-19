#!/usr/bin/env bash
# Run an arbitrary LOCAL Python simulation with installed SCRIMP in Docker.
# Runtime bind mounts DO NOT copy user files into the built image.
# The output directory is local; SCRIMP normally writes ./outputs.
set -euo pipefail

IMAGE="${SCRIMP_DOCKER_IMAGE:-scrimp-mpi:v1.3.0-rc1}"
HOST_DIR="${SCRIMP_SIM_DIR:-$PWD}"
MPI_RANKS="${SCRIMP_MPI_RANKS:-1}"

usage() {
  cat >&2 <<'USAGE'
Usage: bash scripts/run-v130-docker-simulation.sh SCRIPT.py [script arguments...]

Run from the directory containing your local simulations, or set
SCRIMP_SIM_DIR=/absolute/path/to/work-directory.
  SCRIMP_DOCKER_IMAGE=scrimp-mpi:v1.3.0-rc1 (default)
  SCRIMP_MPI_RANKS=1                    (default, or 2/4/etc.)
  SCRIMP_DOCKER_NO_SUDO=1               (only when authorized)
  SCRIMP_MPI_OVERSUBSCRIBE=1            (functional tests only)

The host directory is mounted read-write as /work. SCRIMP's default
./outputs/log, ./outputs/png, ./outputs/mesh, ./outputs/pv remain
in that host directory even after the container is removed.
The container image is not rebuilt or modified.
USAGE
}

if [[ "$#" -lt 1 ]]; then usage; exit 2; fi
if [[ ! -d "$HOST_DIR" ]]; then
  echo "ERROR: simulation working directory does not exist: $HOST_DIR" >&2
  exit 2
fi
HOST_DIR="$(cd "$HOST_DIR" && pwd -P)"
case "$MPI_RANKS" in
  ''|*[!0-9]*) echo "ERROR: SCRIMP_MPI_RANKS must be a positive integer" >&2; exit 2 ;;
esac
if (( MPI_RANKS < 1 )); then
  echo "ERROR: SCRIMP_MPI_RANKS must be >= 1" >&2
  exit 2
fi

if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then
  DOCKER=(docker)
else
  DOCKER=(sudo docker)
fi

# To keep generated files owned by the host user, do not run as the
# image's default mambauser (whose uid/gid can differ). Use the absolute
# conda python / mpiexec binaries to avoid needing shell activation or
# mambauser's home directory with an arbitrary runtime uid.
PYTHON=/opt/conda/envs/scrimp/bin/python
MPIEXEC=/opt/conda/envs/scrimp/bin/mpiexec
RUN_ARGS=(
  run --rm
  --user "$(id -u):$(id -g)"
  --env HOME=/tmp
  --env XDG_CACHE_HOME=/tmp/.cache
  --env PATH=/opt/conda/envs/scrimp/bin:/usr/local/bin:/usr/bin:/bin
  --volume "$HOST_DIR:/work:rw"
  --workdir /work
)

if [[ -n "${DISPLAY:-}" ]]; then
  echo "NOTE: GUI display is not forwarded; plots should be saved to disk."
fi
echo "SCRIMP simulation from $HOST_DIR (mounted /work), image $IMAGE"
echo "Output files are retained in $HOST_DIR/outputs (unless your script overrides SCRIMP paths)."
if (( MPI_RANKS == 1 )); then
  "${DOCKER[@]}" "${RUN_ARGS[@]}" \
    --entrypoint "$PYTHON" \
    "$IMAGE" "$@"
else
  MPI_ARGS=(-n "$MPI_RANKS")
  if [[ "${SCRIMP_MPI_OVERSUBSCRIBE:-0}" == 1 ]]; then
    MPI_ARGS=(--oversubscribe "${MPI_ARGS[@]}")
  fi
  "${DOCKER[@]}" "${RUN_ARGS[@]}" \
    --entrypoint "$MPIEXEC" \
    "$IMAGE" "${MPI_ARGS[@]}" "$PYTHON" "$@"
fi
