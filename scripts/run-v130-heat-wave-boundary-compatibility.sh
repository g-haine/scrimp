#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# Noninvasive original Heat-Wave DAE boundary-compatibility candidates.
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
 echo "Pinned GetFEM/PETSc revision required: $revision" >&2; exit 2;
}
python3 -m py_compile sandbox/probe_v130_heat_wave_boundary_compatibility.py
python3 tests/test_v130_heat_wave_boundary_compatibility.py
python3 tests/test_v130_heat_wave_boundary_compatibility.py
python3 scripts/check-v130-published-immutable.py
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-heat-wave-boundary-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-hw-boundary-candidates-v1","source":"%s","image_revision":"%s"}\n' \
 "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "Original Heat-Wave Gaussian, strict DAE, explicit boundary-DOF candidates"
echo "clamp boundary vs physical mass-minimal L2 zero-trace projection; MPI1"
echo "No change to original source, default initial L2 or DAE tolerance"
echo "No PETSc CN time integration, no physical energy or power PASS"
echo "EVIDENCE: $out"
for candidate in clamp constrained_l2; do
 echo "RUN original fixed Heat-Wave candidate: $candidate"
 if ! "${docker_cmd[@]}" run --rm \
   --user "$(id -u):$(id -g)" --env HOME=/tmp --env MPLBACKEND=Agg \
   --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
   --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
   --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
   --env SCRIMP_MATRIX_MEMORY_MODE=compact \
   --env SCRIMP_GETFEM_RHS_MODE=auto \
   --env SCRIMP_MUMPS_RESET_POLICY=on_change \
   --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
   --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
   --volume "$root/examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro" \
   --volume "$root/sandbox/probe_v130_heat_wave_boundary_compatibility.py:/pilot.py:ro" \
   --volume "$out:/artifacts:rw" --workdir /artifacts \
   "$image" timeout -k 5s "${SCRIMP_V130_BOUNDARY_TIMEOUT:-240s}" \
   mpiexec --oversubscribe --bind-to none -n 1 \
   python -u /pilot.py --candidate "$candidate" \
      --out "/artifacts/$candidate.json" \
   > "$out/$candidate.host.log" 2>&1; then
   echo "FAIL infrastructure executing $candidate; logs preserved: $out" >&2
   tail -n 140 "$out/$candidate.host.log" >&2
   exit 1
 fi
done
python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
root=Path(sys.argv[1])
cases={}
for name in ("clamp","constrained_l2"):
 report=json.loads((root/f"{name}.json").read_text())
 assert report["actual_PETSc_DAE_initializer_executed"] is True
 assert report["no_pseudo_time_or_tolerance_relaxation"] is True
 assert report["T_q_unchanged"] is True
 assert report["after_boundary_Euclidean_coeff_norm"]==0.
 assert report["full_CN_TS_executed"] is False
 assert report["energy_power_qualified"] is False
 assert report["status"] in (
  "PASS_DAE_COMPATIBLE_BOUNDARY_CORRECTION","UNRESOLVED_COMPATIBILITY")
 cases[name]={
  "status":report["status"],
  "corrected_p_relative_coeff_change":
   report["wave_coefficient_relative_change"],
  "native_algebraic_residual":
   report["native_initializer_report"]["native_final_residual_l2"],
  "native_algebraic_tolerance":
   report["native_initializer_report"]["feasibility_tolerance"],
  "mass_optimality":report["field_adjustment_metrics"],
 }
comparison={
 "schema":"scrimp-v130-heat-wave-boundary-candidates-v1",
 "status": (
  "PASS_CONSTRAINED_INITIAL_DAE_FEASIBILITY_ONLY"
  if cases["constrained_l2"]["status"]==
      "PASS_DAE_COMPATIBLE_BOUNDARY_CORRECTION" else
  "UNRESOLVED_CONSTRAINED_INITIAL_DATA"),
 "cases":cases,
 "original_example_modified":False,
 "main_time_integration_executed":False,
 "physical_power_balance_qualified":False,
}
(root/"comparison.json").write_text(json.dumps(comparison,indent=2)+"\n")
with zipfile.ZipFile(Path(str(root)+".zip"),"x",
                     compression=zipfile.ZIP_DEFLATED) as z:
 for path in sorted(root.glob("*")):
  if path.is_file():z.write(path,arcname=path.name)
print(json.dumps(comparison,indent=2))
print("SHAREABLE ORIGINAL HEAT-WAVE BOUNDARY CANDIDATE EVIDENCE:",str(root)+".zip")
if comparison["status"]=="UNRESOLVED_CONSTRAINED_INITIAL_DATA":
 raise SystemExit("Constrained L2 alone did not resolve DAE: inspect native block residuals")
PY
echo "COMPLETE physical boundary compatibility DAE experiment: $out"
