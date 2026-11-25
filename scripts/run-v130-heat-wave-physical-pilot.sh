#!/usr/bin/env bash
# SCRIMP v1.3.0 Heat–Wave 2D: strict native physical energy gate.
# Native MUMPS MPI1 only until DAE + both power balances qualify.
# Always archive native FAIL; never infer physics from PETSc convergence.
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }
docker_cmd=(sudo docker)
[[ "${SCRIMP_DOCKER_NO_SUDO:-0}" != 1 ]] || docker_cmd=(docker)
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
revision="$("${docker_cmd[@]}" image inspect \
  --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == 1d6d0c90115213298047b28237a7f38087563a35 ]] || {
  echo "Wrong pinned GetFEM/PETSc image: $revision" >&2; exit 2;
}
python3 -m py_compile \
  examples/heat_wave.py \
  sandbox/probe_v130_main_2d_examples.py \
  sandbox/measure_v130_heat_wave_physics.py \
  tests/test_v130_main_example_solver_policy.py
python3 tests/test_v130_main_example_solver_policy.py
python3 scripts/check-v130-published-immutable.py
bash -n "$0"
mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-heat-wave-physical-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
mkdir "$out/native" "$out/work"
printf '{"schema":"v130-heat-wave-physical-mumps-mpi1-v1","source":"%s","image_revision":"%s","dt":0.001,"t_final":0.004,"steps":4,"model":"examples/heat_wave.py","MPI":1,"solver":"mumps","full_horizon_qualified":false,"parallel_export_qualified":false}\n' \
  "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "HEAT–WAVE native MUMPS MPI1, four CN steps; physical PASS gated by independent GetFEM balances"
echo "EVIDENCE: $out"
rc=0
"${docker_cmd[@]}" run --rm \
  --user "$(id -u):$(id -g)" \
  --env HOME=/tmp --env MPLBACKEND=Agg \
  --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
  --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
  --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
  --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
  --env SCRIMP_MATRIX_MEMORY_MODE=compact \
  --env SCRIMP_MUMPS_RESET_POLICY=on_change \
  --env SCRIMP_GETFEM_RHS_MODE=auto \
  --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
  --volume "$root/scrimp:/workspace/scrimp/scrimp:ro" \
  --volume "$root/examples/heat_wave.py:/workspace/scrimp/examples/heat_wave.py:ro" \
  --volume "$root/sandbox/probe_v130_main_2d_examples.py:/probe-2d.py:ro" \
  --volume "$root/sandbox/measure_v130_heat_wave_physics.py:/workspace/scrimp/sandbox/measure_v130_heat_wave_physics.py:ro" \
  --volume "$out/native:/artifacts:rw" \
  --volume "$out/work:/scrimp-work:rw" \
  --workdir /scrimp-work \
  "$image" timeout -k 5s "${SCRIMP_V130_HEAT_WAVE_TIMEOUT:-360s}" \
  mpiexec --oversubscribe --bind-to none -n 1 \
  python -u /probe-2d.py --example heat_wave --backend mumps --out /artifacts \
  > "$out/native.host.log" 2>&1 || rc=$?
echo "$rc" > "$out/native.exitcode"
if [[ "$rc" != 0 ]]; then
  echo "UNQUALIFIED native exit=$rc (evidence preserved)" >&2
  tail -n 60 "$out/native.host.log" >&2
fi

python3 - "$out" <<'PY'
import json,sys,zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
import numpy as np

root=Path(sys.argv[1]); native=root/"native"; outputs=root/"work"/"outputs"
rc=int((root/"native.exitcode").read_text())
report={"schema":"v130-heat-wave-physical-audit-v1",
        "status":"FAIL_OR_UNQUALIFIED_HEAT_WAVE_PHYSICAL",
        "source":json.loads((root/"manifest.json").read_text())["source"],
        "native_exit_code":rc,"physical_gates":{},"measurements":{},
        "errors":[],"MPI2_qualified":False,
        "superlu_dist_qualified":False,
        "full_horizon_qualified":False,
        "published_sources_unchanged":True}
try:
    if rc: raise RuntimeError(f"Native MPI1 GetFEM/PETSc exit={rc}")
    meta=json.loads((native/"original_2d_example_report.json").read_text())
    phys=json.loads((native/"heat_wave_physics_diagnostic.json").read_text())
    with np.load(native/"original_2d_example_trace.npz",allow_pickle=False) as z:
        tr={k:z[k].copy() for k in z.files}
    t,z,H=tr["times"],tr["states"],tr["H"]
    if not (meta["example"]=="heat_wave" and meta["factor_solver"]=="mumps"
            and meta["MPI"]==1 and meta["PETSc_reason"]>0
            and meta["steps"]==4 and meta["rejected_steps"]==0
            and meta["SNES_failures"]==0 and meta["heat_wave_fixed_CN"]
            and not meta["heat_wave_physics_certified"]):
        raise RuntimeError("Wrong native fixed-CN diagnostic metadata")
    if not (t.shape==(5,) and H.shape==(5,) and z.ndim==2 and z.shape[0]==5
            and np.all(np.isfinite(H)) and np.all(np.isfinite(z))
            and np.allclose(t,np.arange(5)*.001,atol=1e-11,rtol=0)
            and H[0]>0 and np.max(abs(z-z[0]))>1e-12):
        raise RuntimeError("Invalid real GetFEM FE trajectory or fixed CN times")
    def p(name):
        x=np.asarray(tr["phys_"+name],dtype=float)
        if x.shape!=t.shape or not np.all(np.isfinite(x)):
            raise RuntimeError("Missing/nonfinite physical trace: "+name)
        return x
    scale=max(1.,float(np.max(abs(H))))
    heat,wave,components=p("H_heat"),p("H_wave"),p("H_components")
    d=p("CN_D_heat_flux")[1:]
    interface=p("CN_P_interface_physical")[1:]
    exterior=p("CN_P_outer")[1:]
    dH=np.diff(H); dt=np.diff(t)
    res=dH+dt*(d-interface-exterior)
    he=p("CN_per_step_heat_balance_error")[1:]
    wa=p("CN_per_step_wave_balance_error")[1:]
    total=p("CN_per_step_H_balance_error")[1:]
    cum=p("CN_cumulative_total_balance_error")
    gates={
       "real_fixed_CN_PETSc_steps":True,
       "nonnegative_subsystem_energies":bool(np.min(heat)>=-1e-9*scale and np.min(wave)>=-1e-9*scale),
       "independent_H_component_sum":bool(np.max(abs(H-components))<1e-8*scale),
       "nonnegative_nontrivial_heat_dissipation":bool(np.min(d)>=-1e-10*scale and np.max(d)>1e-10*scale),
       "oriented_interface_power_cancels":bool(np.max(abs(interface))<1e-8*scale),
       "outer_homogeneous_boundary_power_zero":bool(np.max(abs(exterior))<1e-8*scale),
       "CN_heat_balance":bool(np.max(abs(he))<1e-8*scale),
       "CN_wave_balance":bool(np.max(abs(wa))<1e-8*scale),
       "CN_total_independent_balance":bool(max(np.max(abs(res)),np.max(abs(total)),np.max(abs(cum)))<1e-8*scale),
       "closed_H_nonincreasing":bool(np.max(dH)<1e-8*scale),
    }
    report["physical_gates"]=gates
    report["measurements"]={
        "DOFs":int(z.shape[1]),"H_initial":float(H[0]),
        "H_final":float(H[-1]),
        "D_stage_max":float(np.max(d)),
        "max_oriented_interface_power":float(np.max(abs(interface))),
        "max_outer_power":float(np.max(abs(exterior))),
        "max_heat_CN_step_defect":float(np.max(abs(he))),
        "max_wave_CN_step_defect":float(np.max(abs(wa))),
        "max_total_CN_step_defect":float(np.max(abs(res))),
        "max_cumulative_defect":float(np.max(abs(cum))),
        "energy_scale":scale}
    report["diagnostic_source"]=phys
    for var in ("T","p"):
        pvd=outputs/"pv"/var/f"{var}.pvd"
        frames=ET.parse(pvd).findall(".//DataSet")
        if not frames or any(not (pvd.parent/Path(f.attrib["file"])).is_file()
                             for f in frames):
            raise RuntimeError(f"Missing original GetFEM VTK {var} time frames")
    fig=outputs/"png"/"Hamiltonian_Heat1_Wave2_2D.png"
    if not fig.is_file() or fig.stat().st_size<1000:
        raise RuntimeError("Missing original public Hamiltonian PNG")
    if not all(gates.values()):
        raise RuntimeError("Physical gates FAIL: "+",".join(k for k,v in gates.items() if not v))
    report["status"]="PASS_HEAT_WAVE_PHYSICAL_MUMPS_MPI1"
except Exception as exc:
    report["errors"].append(f"{type(exc).__name__}: {exc}")
(root/"summary.json").write_text(json.dumps(report,indent=2)+"\n")
archive=Path(str(root)+".zip")
with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,compresslevel=4) as dest:
    for f in sorted(root.rglob("*")):
        if f.is_file(): dest.write(f,arcname=str(f.relative_to(root)))
print(json.dumps({"status":report["status"],"gates":report["physical_gates"],
                  "measurements":report["measurements"],"errors":report["errors"],
                  "evidence":str(root),"archive":str(archive)},indent=2),flush=True)
raise SystemExit(0 if report["status"].startswith("PASS_") else 1)
PY
