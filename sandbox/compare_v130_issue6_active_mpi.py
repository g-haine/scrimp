# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""MPI #6: independent full original FE/physical/generic-port parity gate.

Original controlled gamma=.25 SWE; the normal feedback does not become
active until t=75. Only t<=.024 tested here. This is a real rotating
wall and nonzero viscous power, not the stationary rest workaround.

Generic Port.compute() returns unweighted CN-stage powers; these are
ONLY compared to identical port quantities on the other MPI ranks.
The original h-weighted SWE powers are checked separately. Never use
plot_powers()'s erroneous trapezoid as a physical balance oracle (#30).
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np

from compare_v130_rotating_ic_parity import compare as checkpoint_compare
from compare_v130_rotating_ic_parity import load_segment


TOL_TIME = 1e-12
TOL_FE_L2 = 1e-8
TOL_OBS = 1e-8
TOL_GENERIC_PORT = 1e-8
PHYSICAL_FIELDS = (
    "H", "volume", "p_normal", "p_viscous",
    "d_grad", "d_div", "boundary_work",
    "viscous_dissipation", "power_balance", "volume_balance",
)
SCHEMA = "scrimp-v130-issue6-active-ports-v1"


def _error(a, b, label, tolerance):
    x, y = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if (x.shape != y.shape or not np.all(np.isfinite(x))
            or not np.all(np.isfinite(y))):
        raise ValueError(f"{label}: nonfinite or mismatched shapes")
    delta = float(np.max(np.abs(x-y)))
    if delta > tolerance:
        raise ValueError(f"{label}: {delta:.9g} > {tolerance:.9g}")
    return delta


def rank_port_traces(directory, ranks, times, fingerprint):
    directory = Path(directory)
    names = {p.name for p in directory.glob("mpi_port_powers_rank*.json")}
    expected = {f"mpi_port_powers_rank{k}.json" for k in range(ranks)}
    if names != expected:
        raise ValueError(
            f"Missing or unexpected rank-local MPI port outputs in {directory}: "
            f"actual={sorted(names)} expected={sorted(expected)}"
        )
    if not (directory/"mpi_physical_powers.png").is_file():
        raise ValueError("Missing rank-0-only h-weighted physical powers plot")
    if (directory/"mpi_physical_powers.png").stat().st_size <= 1000:
        raise ValueError("Truncated rank-0-only h-weighted physical powers plot")
    reports = []
    for rank in range(ranks):
        record = json.loads((
            directory/f"mpi_port_powers_rank{rank}.json"
        ).read_text(encoding="utf-8"))
        if (record.get("schema") != SCHEMA
                or record.get("rank") != rank
                or record.get("world_size") != ranks
                or record.get("source_fingerprint") != fingerprint
                or record.get("physical_h_weighted_balance_is_separate") is not True
                or record.get("normal_feedback_active") is not False
                or "midpoint-stage" not in record.get("sample_semantics", "")):
            raise ValueError(f"Invalid native port metadata rank={rank}")
        _error(record["times"], times, f"rank={rank} physical clock", TOL_TIME)
        powers = record["powers"]
        if not isinstance(powers, dict) or not powers:
            raise ValueError(f"No algebraic port powers on MPI rank {rank}")
        parsed = {}
        for name, values in powers.items():
            vec = np.asarray(values, dtype=float)
            if vec.shape != np.asarray(times).shape or not np.all(np.isfinite(vec)):
                raise ValueError(f"Malformed port trace rank={rank} {name}")
            parsed[name] = vec
        reports.append(parsed)
    names = set(reports[0])
    if any(set(p) != names for p in reports):
        raise ValueError("Port registration differs between MPI ranks")
    intra = 0.
    for rank in range(1,ranks):
        for name in names:
            intra = max(intra,_error(
                reports[rank][name], reports[0][name],
                f"rank={rank} vs rank=0 generic port '{name}'",
                TOL_GENERIC_PORT,
            ))
    activity = max(float(np.max(abs(v)))
                   for p in reports for v in p.values())
    if not math.isfinite(activity) or activity <= 1e-12:
        raise ValueError(
            "Physical rotating SWE MPI port audit is vacuous; "
            "all original generic port powers are zero"
        )
    return reports[0], {"mpi_ranks":ranks, "port_names":sorted(names),
                        "max_intra_MPI_generic_port_difference":intra,
                        "max_original_generic_port_activity":activity,
                        "unique_rank_outputs":len(reports),
                        "root_only_physical_png":True}


def _native_case(directory, phase, ranks):
    report, trace = load_segment(directory, phase)
    if report["MPI"] != ranks or not report.get("mpi_port_power_audit_enabled"):
        raise ValueError("Original MPI #6 native port audit missing or wrong rank count")
    if report["closed_loop_feedback_is_inactive_before_75s"] is not True:
        raise ValueError("Unexpected original feedback clock")
    port, audit=rank_port_traces(
        directory, ranks, trace["times"],
        report["original_FE_fingerprint_sha256"],
    )
    # Original normal/viscous h-weighted power MUST also have nonzero activity.
    if float(np.max(abs(trace["p_viscous"]))) <= 1e-12:
        raise ValueError("Physically rotating wall did not supply nonzero h-weighted power")
    return report, trace, port, audit


def compare_all(seed1, seed4, cont1, cont2, cont4, resume1to4, resume4to1):
    cases = (
        ("seed1",seed1,"seed",1),
        ("seed4",seed4,"seed",4),
        ("continuous1",cont1,"seed",1),
        ("continuous2",cont2,"seed",2),
        ("continuous4",cont4,"seed",4),
        ("resume1to4",resume1to4,"resume",4),
        ("resume4to1",resume4to1,"resume",1),
    )
    evidence={name:_native_case(path,phase,mpi)
              for name,path,phase,mpi in cases}
    ref, a, pref, _ = evidence["continuous1"]
    if not (
        all(x["gamma"]==.25 for x,_,_,_ in evidence.values())
        and len({x["mesh"] for x,_,_,_ in evidence.values()})==1
        and len({x["dt"] for x,_,_,_ in evidence.values()})==1
        and len({x["original_FE_fingerprint_sha256"]
                 for x,_,_,_ in evidence.values()})==1
    ):
        raise ValueError("Original GetFEM physics, FE or absolute time mismatched")
    rows=[]
    for name in ("continuous2","continuous4"):
        report,b,port,_=evidence[name]
        _error(a["times"],b["times"], f"{name} absolute clock",TOL_TIME)
        err_time = np.linalg.norm(a["states"]-b["states"],axis=1)
        err=float(np.max(err_time))
        if not np.isfinite(err) or err > TOL_FE_L2:
            raise ValueError(f"{name} complete FE state MPI parity {err:.9g}")
        observations={
            k:_error(a[k],b[k],f"{name} physical {k}",TOL_OBS)
            for k in PHYSICAL_FIELDS
        }
        if set(port)!=set(pref):
            raise ValueError("Cross MPI generic power port names mismatched")
        generic={
            k:_error(pref[k],port[k],f"{name} generic port {k}",TOL_GENERIC_PORT)
            for k in sorted(port)
        }
        rows.append({
            "case":name, "ranks":report["MPI"],
            "max_complete_FE_trajectory_L2":err,
            "max_original_physical_observable_difference":max(observations.values()),
            "max_original_generic_port_difference":max(generic.values()),
            "per_observable_differences":observations,
            "per_port_differences":generic,
        })
    restart14=checkpoint_compare(seed1,cont4,resume1to4)
    restart41=checkpoint_compare(seed4,cont1,resume4to1)
    return {
        "schema":"scrimp-v130-issue6-compatible-rotating-mpi4-v1",
        "status":"PASS_NATIVE_MPI_1_2_4_ACTIVE_ROTATING_SWE",
        "mesh":ref["mesh"],"gamma":ref["gamma"],
        "dofs":ref["ndof"],"relative_energy_accuracy_not_a_MPI_parity_gate":True,
        "closed_loop_normal_feedback_after_75_not_tested":True,
        "all_original_PETSc_Crank_Nicolson_steps_accepted":True,
        "all_rank_local_port_files_checked":sum(
            audit["unique_rank_outputs"] for _,_,_,audit in evidence.values()
        ),
        "all_7_case_root_only_physical_plots_checked":True,
        "case_port_audits":{name:audit for name,(_,_,_,audit) in evidence.items()},
        "continuous_cross_rank_comparisons":rows,
        "checkpoint_restarts":{"1to4":restart14,"4to1":restart41},
        "max_full_path_FE_L2_MPI":max(r["max_complete_FE_trajectory_L2"] for r in rows),
        "max_generic_port_power_difference_MPI":max(
            r["max_original_generic_port_difference"] for r in rows
        ),
        "max_original_physical_power_observable_difference_MPI":max(
            r["max_original_physical_observable_difference"] for r in rows
        ),
        "max_cold_restart_FE_trajectory_L2":max(
            restart14["maximum_full_FE_trajectory_L2_difference"],
            restart41["maximum_full_FE_trajectory_L2_difference"],
        ),
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("seed1","seed4","cont1","cont2","cont4","resume1to4","resume4to1"):
        p.add_argument("--"+name,required=True,type=Path)
    p.add_argument("--report",required=True,type=Path)
    args=p.parse_args()
    report=compare_all(args.seed1,args.seed4,args.cont1,args.cont2,
                       args.cont4,args.resume1to4,args.resume4to1)
    if args.report.exists():
        raise FileExistsError("Never overwrite original MPI4 native evidence")
    args.report.parent.mkdir(parents=True,exist_ok=True)
    args.report.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(
        f"PASS #6 ACTIVE ROTATING SWE MPI 1/2/4 mesh={report['mesh']} "
        f"DoFs={report['dofs']} "
        f"complete_FE_L2={report['max_full_path_FE_L2_MPI']:.9g} "
        f"rank_local_generic_power_max="
        f"{report['max_generic_port_power_difference_MPI']:.9g} "
        f"physical_power_max="
        f"{report['max_original_physical_power_observable_difference_MPI']:.9g} "
        f"cold_checkpoint_1to4_and_4to1=PASS",
        flush=True,
    )


if __name__=="__main__":
    main()
