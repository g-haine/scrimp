#!/usr/bin/env bash
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.
#
# #30: public Port.compute + DPHS.plot_powers original wave DAE native gate.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
[[ $# == 0 ]] || { echo "Usage: $0" >&2; exit 2; }

python3 -m py_compile \
    scrimp/port.py scrimp/dphs.py scrimp/utils/power_integration.py \
    sandbox/probe_v130_issue30_port_sampling_native.py \
    sandbox/probe_v130_checkpoint_linear_dae.py \
    tests/test_v130_power_integration.py \
    tests/test_v130_port_generic_sampling.py
python3 tests/test_v130_power_integration.py
python3 tests/test_v130_port_generic_sampling.py

docker_exec() {
    if [[ "${SCRIMP_DOCKER_NO_SUDO:-0}" == 1 ]]; then
        docker "$@"
    else
        sudo docker "$@"
    fi
}
image="${SCRIMP_V130_DAE_IMAGE:-scrimp-mpi:v1.3.0-rc3}"
expected="1d6d0c90115213298047b28237a7f38087563a35"
revision="$(docker_exec image inspect --format '{{index .Config.Labels "org.opencontainers.image.revision"}}' "$image")"
[[ "$revision" == "$expected" ]] || {
    echo "Unqualified GetFEM/PETSc image revision=$revision" >&2; exit 2;
}
files=(
    scrimp/dphs.py scrimp/port.py
    scrimp/utils/power_integration.py scrimp/utils/dae_checkpoint.py
    scrimp/utils/mesh.py scrimp/utils/dae_initialization.py
    scrimp/utils/nonlinear_affine_dae.py
    sandbox/probe_v130_issue30_port_sampling_native.py
    sandbox/probe_v130_checkpoint_linear_dae.py
)
for f in "${files[@]}"; do
    [[ -s "$f" ]] || { echo "Missing original native source $f" >&2; exit 2; }
done

mkdir -p sandbox/outputs
out="$(mktemp -d "$root/sandbox/outputs/v130-issue30-port-generic-native-$(git rev-parse --short=8 HEAD)-XXXXXXXX")"
printf '{"schema":"scrimp-v130-issue30-generic-port-native-campaign-v1","source":"%s","image_revision":"%s","original_model":"wave_2D_linear_DAE","schemes":["cn","bdf"],"MPI":[1,2],"native_gate":"PENDING"}\n' \
    "$(git rev-parse HEAD)" "$revision" > "$out/manifest.json"
echo "SCRIMP #30 generic Port.compute ORIGINAL WAVE native smoke"
echo "4 independent original GetFEM/PETSc runs: CN/BDF2 x MPI1/2"
echo "Actual GetFEM port powers, public matplotlib plot, original DAE residual"
echo "CN stage sums vs BDF true instantaneous nodal powers (approx. trapezoid)"
echo "EVIDENCE: $out"

for scheme in cn bdf; do
    for ranks in 1 2; do
        name="${scheme}_mpi${ranks}"
        work="$out/work_$name"
        dest="$out/$name"
        log="$out/$name.host.log"
        mkdir -p "$work" "$dest"
        echo "RUN scheme=$scheme MPI=$ranks original physical wave"
        if ! docker_exec run --rm \
            --user "$(id -u):$(id -g)" \
            --workdir /scrimp-work \
            --env HOME=/tmp --env MPLBACKEND=Agg \
            --env PYTHONPATH=/workspace/scrimp:/workspace/scrimp/sandbox \
            --env OMP_NUM_THREADS=1 --env OPENBLAS_NUM_THREADS=1 \
            --env SCRIMP_GETFEM_MATRIX_MODE=distributed \
            --env SCRIMP_GETFEM_BRIDGE_DISTRIBUTED=coo \
            --env SCRIMP_MATRIX_MEMORY_MODE=compact \
            --env SCRIMP_MUMPS_RESET_POLICY=on_change \
            --env SCRIMP_GETFEM_RHS_MODE=auto \
            --env SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem \
            --volume "$work:/scrimp-work:rw" \
            --volume "$dest:/artifacts:rw" \
            --volume "$root/scrimp/dphs.py:/workspace/scrimp/scrimp/dphs.py:ro" \
            --volume "$root/scrimp/port.py:/workspace/scrimp/scrimp/port.py:ro" \
            --volume "$root/scrimp/utils/power_integration.py:/workspace/scrimp/scrimp/utils/power_integration.py:ro" \
            --volume "$root/scrimp/utils/dae_checkpoint.py:/workspace/scrimp/scrimp/utils/dae_checkpoint.py:ro" \
            --volume "$root/scrimp/utils/mesh.py:/workspace/scrimp/scrimp/utils/mesh.py:ro" \
            --volume "$root/scrimp/utils/dae_initialization.py:/workspace/scrimp/scrimp/utils/dae_initialization.py:ro" \
            --volume "$root/scrimp/utils/nonlinear_affine_dae.py:/workspace/scrimp/scrimp/utils/nonlinear_affine_dae.py:ro" \
            --volume "$root/sandbox/probe_v130_checkpoint_linear_dae.py:/workspace/scrimp/sandbox/probe_v130_checkpoint_linear_dae.py:ro" \
            --volume "$root/sandbox/probe_v130_issue30_port_sampling_native.py:/workspace/scrimp/sandbox/probe_v130_issue30_port_sampling_native.py:ro" \
            "$image" timeout -k 5s "${SCRIMP_V130_ISSUE30_TIMEOUT:-360s}" \
            mpiexec --oversubscribe --bind-to none -n "$ranks" \
            python -u /workspace/scrimp/sandbox/probe_v130_issue30_port_sampling_native.py \
                --scheme "$scheme" --out /artifacts \
            > "$log" 2>&1; then
            echo "FAIL original native public plot $name: evidence retained $out" >&2
            tail -n 140 "$log" >&2
            exit 1
        fi
        grep -a '^NATIVE PUBLIC PLOT POWER PASS' "$log" || {
            echo "FAIL missing verified native completion $name" >&2
            tail -n 100 "$log" >&2
            exit 1
        }
        [[ -s "$dest/native_public_power_report.json" \
            && -s "$dest/original_physical_wave_power.npz" \
            && -s "$dest/public_plot_powers.png" ]] || {
            echo "FAIL incomplete original public-plot evidence $name" >&2
            exit 1
        }
    done
done

python3 - "$out" <<'PY'
import json
import sys
import zipfile
from pathlib import Path
import numpy as np

out=Path(sys.argv[1]).resolve()
all_cases={}
for scheme in ("cn","bdf"):
    for mpi in (1,2):
        key=f"{scheme}_mpi{mpi}"
        directory=out/key
        report=json.loads((directory/"native_public_power_report.json").read_text())
        expected="cn_step" if scheme=="cn" else "instantaneous"
        assert report["status"]=="NATIVE_PUBLIC_PORT_COMPUTE_AND_PLOT_POWERS_PASS"
        assert report["scheme"]==scheme and report["sampling"]==expected
        assert report["mpi_ranks"]==mpi and report["steps"]==12
        assert report["getfem_algebraic_residual_max_l2"]<1e-8
        assert report["maximum_rank_local_port_difference"]<1e-8
        assert report["maximum_original_generic_port_activity"]>1e-9
        if scheme=="cn":
            assert report["maximum_CN_per_step_signed_energy_defect"]<1e-11
        else:
            assert report["BDF_nodal_plot_is_approximate_not_discrete_quadrature"] is True
        assert (directory/"public_plot_powers.png").stat().st_size>1000
        with np.load(directory/"original_physical_wave_power.npz",allow_pickle=False) as a:
            data={name:a[name] for name in a.files}
        assert data["times"].shape==(13,)
        assert data["states"].shape[0]==13
        assert data["raw_ports"].shape[1]==13
        assert data["signed_port_integrals"].shape==data["raw_ports"].shape
        assert all(np.all(np.isfinite(x)) for x in data.values() if x.dtype.kind in "fiu")
        all_cases[key]=(report,data)
results={}
for scheme in ("cn","bdf"):
    x=all_cases[f"{scheme}_mpi1"][1]
    y=all_cases[f"{scheme}_mpi2"][1]
    assert np.array_equal(x["port_names"],y["port_names"])
    err=float(np.max(np.linalg.norm(x["states"]-y["states"],axis=1)))
    diff=float(np.max(abs(x["raw_ports"]-y["raw_ports"])))
    work=float(np.max(abs(x["signed_port_integrals"]-y["signed_port_integrals"])))
    assert err<1e-8 and diff<1e-8 and work<1e-8, (scheme,err,diff,work)
    results[scheme]={"MPI1_vs_MPI2_max_full_FE_trajectory_L2":err,
                     "MPI1_vs_MPI2_max_raw_port_power_difference":diff,
                     "MPI1_vs_MPI2_max_public_integrated_port_difference":work}
summary={
    "schema":"scrimp-v130-issue30-port-generic-native-campaign-v1",
    "native_gate":"PASS_PUBLIC_PORT_COMPUTE_GENERIC_CN_AND_BDF_PLOTTING",
    "source":json.loads((out/"manifest.json").read_text())["source"],
    "independent_original_wave_PETSc_runs":4,
    "native_CN_and_BDF_Mpi1_2_cases":list(all_cases),
    "results":results,
    "BDF_does_not_claim_discrete_stage_quadrature":True,
    "original_GetFEM_weak_form_and_PETSc_solves_unchanged":True,
}
(out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2))
zip_path=Path(str(out)+".zip")
with zipfile.ZipFile(zip_path,"x",compression=zipfile.ZIP_DEFLATED) as z:
    for p in sorted(out.rglob("*")):
        if p.is_file() and ("work_" not in str(p.relative_to(out))):
            z.write(p,arcname=str(p.relative_to(out)))
print("SHAREABLE ORIGINAL NATIVE EVIDENCE ZIP:",zip_path)
PY
echo "COMPLETE original native GetFEM/PETSc port-generic evidence: $out"
