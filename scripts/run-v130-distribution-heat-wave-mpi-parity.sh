#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

# Bounded real Heat-Wave PDE MPI1/MPI2 parity INSIDE the released image.
# No SCRIMP source-tree mount. Native FE Hamiltonians, not generated data.
# This does not certify MPI PDEs on every mesh, topology or time horizon.
set -Eeuo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

image="${SCRIMP_V130_DISTRIBUTION_IMAGE:-}"
expected="${SCRIMP_V130_EXPECT_REVISION:-}"
if [[ -z "$image" || "$image" == *:latest || "$image" == "scrimp-mpi:local" ]]; then
  echo "ERROR: set SCRIMP_V130_DISTRIBUTION_IMAGE to a qualified candidate tag" >&2
  exit 2
fi
if [[ ! "$expected" =~ ^[0-9a-f]{40}$ ]]; then
  echo "ERROR: set SCRIMP_V130_EXPECT_REVISION to the image's exact Git SHA" >&2
  exit 2
fi
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then
  docker_cmd=(docker)
else
  docker_cmd=(sudo docker)
fi

revision="$("${docker_cmd[@]}" image inspect --format \
    '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
platform="$("${docker_cmd[@]}" image inspect --format '{{.Os}}/{{.Architecture}}' "$image")"
if [[ "$revision" != "$expected" ]]; then
  echo "FAIL: image revision $revision differs from requested $expected" >&2
  exit 3
fi
if [[ "$platform" != "linux/amd64" ]]; then
  echo "FAIL: this qualification is scoped to linux/amd64, got $platform" >&2
  exit 3
fi

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-distribution-heat-wave-mpi-XXXXXXXX")"
exec > >(tee "$out/gate.log") 2>&1
echo "SCRIMP v1.3.0 native Heat-Wave PDE MPI1/MPI2 distribution gate"
echo "Image: $image; revision: $revision; platform: $platform"
echo "Evidence: $out"
echo "Only the independent test/evidence directory is mounted at runtime."
echo "No library, examples, or SCRIMP repository source mount."

cat > "$out/probe_heat_wave.py" <<'PY'
"""Execute a real canonical coupled Heat-Wave PDE from the installed image."""
import json
import os
import warnings
from pathlib import Path

import numpy as np
from petsc4py import PETSc
from examples.heat_wave import heat_wave_eq

warnings.simplefilter("error", np.exceptions.ComplexWarning)
comm = PETSc.COMM_WORLD
rank = comm.getRank()
size = comm.getSize()
expected_ranks = int(os.environ["SCRIMP_MPI_PROBE_N"])
assert size == expected_ranks, (size, expected_ranks)

# Identical native weak forms, mesh, constrained initial data and MUMPS LU.
# physics_only suppresses rank-unsafe plots/VTK writes but NO physics.
model = heat_wave_eq(
    final_time=0.02, time_step=0.002, save_step=0.002,
    initial_wave_boundary="clamp", initial_wave_profile="gaussian",
    physics_only=True,
)
times = np.asarray(model.solution["t"], dtype=float).ravel()
energy = np.asarray(model.get_Hamiltonian(), dtype=float).ravel()
assert times.shape == energy.shape == (11,), (times.shape, energy.shape)
assert np.isfinite(times).all() and np.isfinite(energy).all()
assert np.all(np.diff(times) > 0), times
assert abs(times[0]) < 1e-10 and abs(times[-1] - 0.02) < 1e-10
assert np.all(energy >= -1e-10), energy
assert energy[0] > 0.0
comm.barrier()
if rank == 0:
    target = Path.cwd() / "series.json"
    target.write_text(json.dumps({
        "schema": "scrimp-v130-native-heat-wave-mpi-energy-v1",
        "ranks": size, "method": "original GetFEM Hamiltonian terms",
        "time": times.tolist(), "Hamiltonian": energy.tolist(),
        "pde": "examples.heat_wave.heat_wave_eq",
        "physics_only": True, "no_field_interpolation": True,
    }, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: native {size}-rank Heat-Wave PDE, "
          f"{len(times)} genuine FE energy samples -> {target}", flush=True)
comm.barrier()
PY

for ranks in 1 2; do
  mkdir -p "$out/mpi$ranks"
  echo "=== Running original Heat-Wave PDE with $ranks MPI rank(s) ==="
  "${docker_cmd[@]}" run --rm --user "$(id -u):$(id -g)" \
    --env HOME=/tmp --env XDG_CACHE_HOME=/tmp/.cache \
    --env MPLBACKEND=Agg --env OMP_NUM_THREADS=1 \
    --env OPENBLAS_NUM_THREADS=1 --env PYTHONPATH=/workspace/scrimp \
    --env SCRIMP_MPI_PROBE_N="$ranks" \
    --mount "type=bind,src=$out,dst=/work" \
    --workdir "/work/mpi$ranks" \
    --entrypoint /opt/conda/envs/scrimp/bin/mpiexec \
    "$image" --oversubscribe --bind-to none -n "$ranks" \
    /opt/conda/envs/scrimp/bin/python -u /work/probe_heat_wave.py \
    2>&1 | tee "$out/mpi$ranks/run.log"
  test -s "$out/mpi$ranks/series.json" || {
    echo "FAIL: no genuine native FE trajectory from MPI=$ranks" >&2
    exit 3
  }
done

python3 - "$out" "$image" "$revision" "$platform" <<'PY'
import json
import math
import sys
from pathlib import Path

root = Path(sys.argv[1])
one = json.loads((root / "mpi1/series.json").read_text(encoding="utf-8"))
two = json.loads((root / "mpi2/series.json").read_text(encoding="utf-8"))
assert (one["ranks"], two["ranks"]) == (1, 2)
t1, t2 = one["time"], two["time"]
h1, h2 = one["Hamiltonian"], two["Hamiltonian"]
assert len(t1) == len(t2) == len(h1) == len(h2) == 11
assert all(math.isfinite(x) for x in t1+t2+h1+h2)
assert max(abs(a-b) for a,b in zip(t1,t2)) < 1e-10
norm_ref = math.sqrt(sum(x*x for x in h1))
assert norm_ref > 0
relative_l2 = math.sqrt(sum((a-b)**2 for a,b in zip(h1,h2))) / norm_ref
tolerance = 1e-7
if relative_l2 > tolerance:
    raise AssertionError(
        f"Heat-Wave native MPI1/MPI2 physical H parity failed: "
        f"rel_l2={relative_l2:.9g} > {tolerance}"
    )
report = {
    "schema": "scrimp-v130-distribution-heat-wave-mpi-gate-v1",
    "status": "PASS_NATIVE_HEAT_WAVE_MPI1_MPI2_HAMILTONIAN_PARITY",
    "image": sys.argv[2], "source_revision": sys.argv[3],
    "platform": sys.argv[4],
    "native_FE_saved_samples_per_case": 11,
    "time_final": 0.02, "dt": 0.002,
    "relative_l2_H_mpi1_vs_mpi2": relative_l2,
    "parity_tolerance": tolerance,
    "source_isolation": "image-internal SCRIMP and canonical example; test-only host mount",
    "scope": "one real original Heat-Wave PDE, bounded native energy trajectory parity",
    "not_tested": [
        "pointwise or full-state FE MPI equivalence",
        "large/long PDE runs or MPI>2",
        "fresh-machine image pull",
        "cross-platform portability",
    ],
}
(root / "summary.json").write_text(
    json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
print(json.dumps(report, indent=2, sort_keys=True))
PY
echo "PASS: native Heat-Wave distribution PDE MPI1/MPI2 energy parity."
