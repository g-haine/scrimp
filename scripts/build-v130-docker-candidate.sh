#!/usr/bin/env bash
# Build SCRIMP v1.3.0 candidate under a NEW tag, never scrimp-mpi:local.
# BuildKit retains apt, micromamba, pip and official GetFEM source caches.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${SCRIMP_DOCKER_CANDIDATE_IMAGE:-scrimp-mpi:v1.3.0-rc1}"
MAKE_JOBS="${MAKE_JOBS:-4}"
# The release candidate is built from GitHub by default, not from a
# simulation-filled developer checkout. "local" is an explicit dev opt-in.
SCRIMP_DOCKER_SOURCE="${SCRIMP_DOCKER_SOURCE:-github}"
SCRIMP_GITHUB_REPO="${SCRIMP_GITHUB_REPO:-https://github.com/g-haine/scrimp.git}"
SCRIMP_GITHUB_REF="${SCRIMP_GITHUB_REF:-release/1.3.0-preparation}"
case "$SCRIMP_DOCKER_SOURCE" in
  github|local) ;;
  *) echo "ERROR: SCRIMP_DOCKER_SOURCE must be github or local" >&2; exit 2 ;;
esac
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

echo "SCRIMP v1.3.0 candidate BuildKit build" | tee "$LOG"
echo "Image: $IMAGE; mode=$SCRIMP_DOCKER_SOURCE; MAKE_JOBS=$MAKE_JOBS" | tee -a "$LOG"
echo "Qualified image scrimp-mpi:local remains untouched." | tee -a "$LOG"

# The LIVE repository must NEVER be the Docker context. The old
# .dockerignore negation rules accidentally transferred >45 GB. Even
# a GitHub checkout is filtered into a physically distinct tiny context.
STAGED_CONTEXT="$(mktemp -d -t scrimp-v130-context.XXXXXXXX)"
SOURCE_TEMP=""
cleanup() {
  rm -rf -- "$STAGED_CONTEXT"
  if [[ -n "$SOURCE_TEMP" ]]; then rm -rf -- "$SOURCE_TEMP"; fi
}
trap cleanup EXIT

if [[ "$SCRIMP_DOCKER_SOURCE" == github ]]; then
  # Fetch one immutable Git object graph snapshot from GitHub. This
  # starts from an EMPTY /tmp directory, independent of host simulations.
  SOURCE_TEMP="$(mktemp -d -t scrimp-v130-github.XXXXXXXX)"
  echo "Fetching GitHub source $SCRIMP_GITHUB_REPO @ $SCRIMP_GITHUB_REF" \
    | tee -a "$LOG"
  git init -q "$SOURCE_TEMP"
  git -C "$SOURCE_TEMP" remote add origin "$SCRIMP_GITHUB_REPO"
  GIT_LFS_SKIP_SMUDGE=1 git -C "$SOURCE_TEMP" fetch --quiet \
    --no-tags --depth 1 origin "$SCRIMP_GITHUB_REF"
  GIT_LFS_SKIP_SMUDGE=1 git -C "$SOURCE_TEMP" checkout --detach --quiet FETCH_HEAD
  SOURCE_TREE="$SOURCE_TEMP"
else
  echo "Local Git source explicitly selected (development only)" | tee -a "$LOG"
  SOURCE_TREE="$ROOT"
  if [[ -n "$(git -C "$SOURCE_TREE" status --porcelain --untracked-files=no)" ]]; then
    echo "WARNING: local tracked changes are included, not represented by the Git SHA." \
      | tee -a "$LOG"
  fi
fi

REVISION="$(git -C "$SOURCE_TREE" rev-parse HEAD)"
echo "Exact Git revision: $REVISION" | tee -a "$LOG"
python3 "$SOURCE_TREE/scripts/audit-v130-docker-context.py" \
  --stage-to "$STAGED_CONTEXT" 2>&1 | tee -a "$LOG"
echo "Build context: physically isolated $STAGED_CONTEXT, from $SCRIMP_DOCKER_SOURCE" \
  | tee -a "$LOG"

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
  --file "$STAGED_CONTEXT/docker/Dockerfile" \
  --tag "$IMAGE" \
  --build-arg "SCRIMP_REVISION=$REVISION" \
  --build-arg "MAKE_JOBS=$MAKE_JOBS" \
  "$STAGED_CONTEXT" 2>&1 | tee -a "$LOG"

echo "Built image ID and revision:" | tee -a "$LOG"
"${DOCKER[@]}" image inspect --format \
  '{{.Id}}  source-revision={{index .Config.Labels "org.opencontainers.image.revision"}}' \
  "$IMAGE" | tee -a "$LOG"

echo "NEXT: bash scripts/smoke-v130-docker-candidate.sh" | tee -a "$LOG"
