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

# Use Buildx when available (including Docker Desktop); otherwise use the
# standard BuildKit builder. No implicit push, pruning or Docker image deletion.
if "${DOCKER[@]}" buildx version >/dev/null 2>&1; then
  "${DOCKER[@]}" buildx build --load \
    --progress=plain \
    --file docker/Dockerfile \
    --tag "$IMAGE" \
    --build-arg "SCRIMP_REVISION=$REVISION" \
    --build-arg "MAKE_JOBS=$MAKE_JOBS" \
    . 2>&1 | tee -a "$LOG"
else
  # A legacy Docker daemon may not enable BuildKit by default.
  # Force BuildKit because the Dockerfile uses RUN cache mounts.
  if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then
    DOCKER_BUILDKIT=1 docker build --progress=plain \
      --file docker/Dockerfile \
      --tag "$IMAGE" \
      --build-arg "SCRIMP_REVISION=$REVISION" \
      --build-arg "MAKE_JOBS=$MAKE_JOBS" \
      . 2>&1 | tee -a "$LOG"
  else
    sudo env DOCKER_BUILDKIT=1 docker build --progress=plain \
      --file docker/Dockerfile \
      --tag "$IMAGE" \
      --build-arg "SCRIMP_REVISION=$REVISION" \
      --build-arg "MAKE_JOBS=$MAKE_JOBS" \
      . 2>&1 | tee -a "$LOG"
  fi
fi

echo "Built image ID and revision:" | tee -a "$LOG"
"${DOCKER[@]}" image inspect --format \
  '{{.Id}}  source-revision={{index .Config.Labels "org.opencontainers.image.revision"}}' \
  "$IMAGE" | tee -a "$LOG"

echo "NEXT: bash scripts/smoke-v130-docker-candidate.sh" | tee -a "$LOG"
