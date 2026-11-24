#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Original Heat-Wave DAE feasibility REJECTION diagnostic (not relaxation).
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }
docker_cmd=(sudo docker)
if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then docker_cmd=(docker); fi
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
  --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "1d6d0c90115213298047b28237a7f38087563a35" ]] || {
 echo "Unqualified GetFEM/PETSc image: $revision" >&2; exit 2;
}
python3 -m py_compile sandbox/probe_v130_heat_wave_dae_feasibility.py
python3 tests/test_v130_dae_algebraic_block_diagnostics.py
python3 tests/test_v130_dae_algebraic_block_diagnostics.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-heat-wave-dae-feasibility-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-original-heat-wave-dae-feasibility-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "Original Heat-Wave fixed differential Gaussian, default L2, COO, MPI1"
echo "Expect correct DAE inconsistency rejection; record actual residual by algebraic port"
echo "No pseudo time, no altered model, no changed DAE thresholds"
echo "EVIDENCE: $out"
if ! "${docker_cmd[@]}" run --rm \
 --user "$(id -u):$(id -g)" \
 --env HOME=/tmp --env MPLBACKEND=Agg \
 --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
 --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
 --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
 --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
 --env SCRIMP_GETFEM_RHS_MODE=auto \
 --env SCRIMP_MATRIX_MEMORY_MODE=compact \
 --env SCRIMP_MUMPS_RESET_POLICY=on_change \
 --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
 --volume "$root/examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro" \
 --volume "$root/sandbox/probe_v130_heat_wave_dae_feasibility.py:/workspace/scrimp/sandbox/probe_v130_heat_wave_dae_feasibility.py:ro" \
 --volume "$out:/artifacts:rw" \
 --workdir /artifacts \
 "$image" timeout -k 5s "${SCRIMP_V130_DAE_FEASIBILITY_TIMEOUT:-240s}" \
 mpiexec --oversubscribe --bind-to none -n 1 \
 python -u /workspace/scrimp/sandbox/probe_v130_heat_wave_dae_feasibility.py \
 --out /artifacts/feasibility.json \
 > "$out/native.host.log" 2>&1; then
 echo "FAIL expected-incompatibility diagnostics, original logs: $out" >&2
 tail -n 150 "$out/native.host.log" >&2; exit 1
fi
grep -a '^PASS_DETECTED_PHYSICALLY_INCOMPATIBLE_INITIAL_DATA$' "$out/native.host.log"
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1])
r=json.loads((root/"feasibility.json").read_text())
assert r["status"]=="PASS_DETECTED_PHYSICALLY_INCOMPATIBLE_INITIAL_DATA"
assert r["full_TS_executed"] is False
assert r["physical_model_qualified"] is False
assert r["native_initializer_report"]["status"]=="UNRESOLVED"
assert r["native_initializer_report"]["differential_initial_exact"] is True
assert r["analytical_wave_gaussian_at_outer_boundary_x1_y0"]>0.09
print(json.dumps(r,indent=2))
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED) as z:
 for f in sorted(root.glob("*")):
  if f.is_file(): z.write(f,arcname=f.name)
print("SHAREABLE ORIGINAL HEAT-WAVE FEASIBILITY DIAGNOSTICS:",archive)
PY
echo "COMPLETE original Heat-Wave consistent-data diagnostic: $out"
