#!/usr/bin/env bash
# Build SCRIMP v1.3.0 candidate under a NEW tag, never scrimp-mpi:local.
# BuildKit retains apt, micromamba, pip and official GetFEM source caches.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${SCRIMP_DOCKER_CANDIDATE_IMAGE:-scrimp-mpi:v1.3.0-rc1}"
MAKE_JOBS="${MAKE_JOBS:-4}"
case "$IMAGE" in
  scrimp-mpi:local|scrimp-mpi:latest|"")
    echo "ERROR: refusing to overwrite the qualified scrimp-mpi:local image." >&2
    exit 2
    ;;
esac

if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then
  DOCKER=(docker)
else
  DOCKER=(sudo docker)
fi

cd "$ROOT"
mkdir -p sandbox/outputs
LOG="$ROOT/sandbox/outputs/v130-docker-build.log"
REVISION="$(git rev-parse HEAD 2>/dev/null || echo unknown)"

echo "SCRIMP v1.3.0 candidate BuildKit build" | tee "$LOG"
echo "Image: $IMAGE; source: $REVISION; MAKE_JOBS=$MAKE_JOBS" | tee -a "$LOG"
echo "Qualified image scrimp-mpi:local remains untouched." | tee -a "$LOG"

# A BuildKit Dockerfile needs the Buildx CLI plugin. In particular,
# DOCKER_BUILDKIT=1 is NOT a valid fallback when Buildx is missing: recent
# Docker versions abort with "buildx component is missing or broken".
# Probe through the same identity (sudo or non-sudo) used for the build.
if ! "${DOCKER[@]}" buildx version 2>&1 | tee -a "$LOG"; then
  {
    echo "ERROR: Docker Buildx is missing or not accessible through ${DOCKER[*]}."
    echo "Ubuntu docker.io installation: sudo apt-get install docker-buildx"
    echo "Official Docker CE apt installation: sudo apt-get install docker-buildx-plugin"
    echo "First inspect: apt-cache policy docker.io docker-ce docker-buildx docker-buildx-plugin"
    echo "Then verify: ${DOCKER[*]} buildx version"
    echo "This script will NOT install packages or replace Docker automatically."
  } | tee -a "$LOG" >&2
  exit 2
fi

# Buildx is explicitly required. Never push, prune or delete images.
"${DOCKER[@]}" buildx build --load \
  --progress=plain \
  --file docker/Dockerfile \
  --tag "$IMAGE" \
  --build-arg "SCRIMP_REVISION=$REVISION" \
  --build-arg "MAKE_JOBS=$MAKE_JOBS" \
  . 2>&1 | tee -a "$LOG"

echo "Built image ID and revision:" | tee -a "$LOG"
"${DOCKER[@]}" image inspect --format \
  '{{.Id}}  source-revision={{index .Config.Labels "org.opencontainers.image.revision"}}' \
  "$IMAGE" | tee -a "$LOG"

echo "NEXT: bash scripts/smoke-v130-docker-candidate.sh" | tee -a "$LOG"
