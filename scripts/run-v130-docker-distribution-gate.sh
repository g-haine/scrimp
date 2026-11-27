#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# SCRIMP v1.3.0: fail-closed, source-independent USER runtime qualification.
# Build/packaging happen separately; this never pulls/pushes/prunes images.
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${SCRIMP_V130_DISTRIBUTION_IMAGE:-}"
if [[ -z "$IMAGE" || "$IMAGE" == *:latest || "$IMAGE" == scrimp-mpi:local ]]; then
    echo "ERROR: set SCRIMP_V130_DISTRIBUTION_IMAGE to a NEW, explicit candidate tag." >&2
    exit 2
fi
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then
    DOCKER=(docker)
else
    DOCKER=(sudo docker)
fi

# The candidate must have been built from a known revision. Do not accept
# a mutable branch name or silently use the locally cached image ID alone.
EXPECTED="${SCRIMP_V130_EXPECT_REVISION:-}"
if [[ ! "$EXPECTED" =~ ^[0-9a-f]{40}$ ]]; then
    echo "ERROR: set SCRIMP_V130_EXPECT_REVISION to the exact 40-char Git SHA." >&2
    exit 2
fi
# The release ZIP must be created from exactly the code version that
# produced the candidate image; do not mix two Git revisions.
SOURCE_SHA="$(git -C "$ROOT" rev-parse HEAD)"
if [[ "$SOURCE_SHA" != "$EXPECTED" ]]; then
    echo "ERROR: workspace revision $SOURCE_SHA != candidate source $EXPECTED" >&2
    exit 3
fi
if [[ -n "$(git -C "$ROOT" status --porcelain --untracked-files=no)" ]]; then
    echo "ERROR: tracked local modifications would make the ZIP untraceable." >&2
    exit 3
fi
mkdir -p "$ROOT/sandbox/outputs"
OUT="$(mktemp -d "$ROOT/sandbox/outputs/v130-docker-distribution-XXXXXXXX")"
exec > >(tee "$OUT/gate.log") 2>&1
echo "SCRIMP v1.3.0 Docker distribution gate"
echo "Evidence: $OUT"
echo "Image: $IMAGE; expected source: $EXPECTED"
echo "No image is pulled, published, modified or deleted."

IMAGE_ID="$("${DOCKER[@]}" image inspect --format '{{.Id}}' "$IMAGE")"
REVISION="$("${DOCKER[@]}" image inspect --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$IMAGE")"
PLATFORM="$("${DOCKER[@]}" image inspect --format '{{.Os}}/{{.Architecture}}' "$IMAGE")"
if [[ "$REVISION" != "$EXPECTED" ]]; then
    echo "FAIL: candidate source label ($REVISION) != expected SHA ($EXPECTED)" >&2
    exit 3
fi
echo "Image ID: $IMAGE_ID; platform: $PLATFORM; revision: $REVISION"

# Metadata must agree with docs/source/conf.py and setup.py on this branch.
# Do not mistake an older image containing an editable 1.2.0 package for
# the newly prepared distribution image.
"${DOCKER[@]}" run --rm --entrypoint /opt/conda/envs/scrimp/bin/python "$IMAGE" -c '
from importlib.metadata import version
actual = version("scrimp")
assert actual == "1.3.0", ("wrong wheel/install version", actual)
print("PASS: installed SCRIMP distribution version", actual)
'

# Use the same ZIP packaging machinery as GitHub Releases, then extract it
# into an isolated user folder. No source-tree volume or SCRIMP checkout is
# ever mounted into the candidate for the rest of this test.
python3 "$ROOT/scripts/package-scrimp-project-template.py" \
    --output "$OUT/scrimp-project-template.zip" --image "$IMAGE"
python3 - "$OUT/scrimp-project-template.zip" "$OUT" "$IMAGE" <<'PY'
from pathlib import Path
from zipfile import ZipFile
import json
import sys

archive, out, image = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
with ZipFile(archive) as z:
    assert z.testzip() is None, "invalid ZIP CRC"
    members = z.namelist()
    assert members and all(
        p.startswith("scrimp-project/") and
        ".." not in Path(p).parts and
        not Path(p).is_absolute()
        for p in members
    ), "unsafe or unexpected ZIP member"
    assert "scrimp-project/simulations/first_simulation.py" in members
    assert "scrimp-project/simulations/mpi_check.py" in members
    assert not any("sandbox/" in p or ".git/" in p for p in members)
    config = json.loads(z.read("scrimp-project/.devcontainer/devcontainer.json"))
    assert config["image"] == image, (config["image"], image)
    z.extractall(out / "user-workspace")
assert (out / "user-workspace/scrimp-project/outputs").is_dir()
print("PASS: standalone user ZIP, CRC, pinned image and extraction")
PY
PROJECT="$OUT/user-workspace/scrimp-project"
echo "=== Independent user simulation: real host-persistent output ==="
"${DOCKER[@]}" run --rm --user "$(id -u):$(id -g)" \
    --env HOME=/tmp --env XDG_CACHE_HOME=/tmp/.cache \
    --mount "type=bind,src=$PROJECT,dst=/work" --workdir /work \
    --entrypoint /opt/conda/envs/scrimp/bin/python "$IMAGE" \
    simulations/first_simulation.py

python3 - "$PROJECT/outputs/oscillator.csv" <<'PY'
import csv
import math
import sys
from pathlib import Path

csv_file = Path(sys.argv[1])
assert csv_file.is_file() and csv_file.stat().st_size > 0, csv_file
with csv_file.open(newline="", encoding="utf-8") as stream:
    reader = csv.DictReader(stream)
    assert reader.fieldnames == ["time", "q", "p", "H"], reader.fieldnames
    rows = list(reader)
assert len(rows) == 11, len(rows)
times = [float(row["time"]) for row in rows]
assert abs(times[0]) < 1e-12 and abs(times[-1] - 0.1) < 1e-10, times
assert all(b > a for a, b in zip(times, times[1:])), times
for row in rows:
    q, p, energy = [float(row[x]) for x in ("q", "p", "H")]
    assert all(math.isfinite(x) for x in (q, p, energy)), row
    assert math.isclose(energy, 2.5*q*q + p*p/3.0,
                        rel_tol=1e-10, abs_tol=1e-12), row
print("PASS: user-owned, persistent host CSV and real oscillator energy")
PY

echo "=== Independent user workspace: two MPI ranks ==="
"${DOCKER[@]}" run --rm --user "$(id -u):$(id -g)" \
    --env HOME=/tmp --env XDG_CACHE_HOME=/tmp/.cache \
    --mount "type=bind,src=$PROJECT,dst=/work" --workdir /work \
    --entrypoint /opt/conda/envs/scrimp/bin/mpiexec "$IMAGE" \
    --oversubscribe -n 2 /opt/conda/envs/scrimp/bin/python simulations/mpi_check.py

python3 - "$OUT" "$IMAGE" "$IMAGE_ID" "$REVISION" "$PLATFORM" <<'PY'
from pathlib import Path
import json
import sys

out = Path(sys.argv[1])
report = {
    "schema": "scrimp-v130-docker-distribution-gate-v1",
    "status": "PASS_STANDALONE_ZIP_IMAGE_SERIAL_MPI_STARTUP",
    "image": sys.argv[2],
    "image_id": sys.argv[3],
    "source_revision": sys.argv[4],
    "platform": sys.argv[5],
    "distribution_version": "1.3.0",
    "zip": str(out / "scrimp-project-template.zip"),
    "user_csv": str(out / "user-workspace/scrimp-project/outputs/oscillator.csv"),
    "scope": "native metadata, standalone ZIP, serial physics + host output, MPI startup only",
    "not_tested": [
        "fresh host or foreign-machine Docker daemon",
        "numerical PDE MPI physical parity",
        "VS Code GUI F5",
        "Windows/macOS",
        "GHCR publishing or registry pulls",
    ],
}
(out / "summary.json").write_text(
    json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps(report, indent=2))
PY
echo "PASS: bounded Docker distribution gate, no image publication."
