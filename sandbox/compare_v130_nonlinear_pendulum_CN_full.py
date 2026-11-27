"""Independent native evidence check for canonical PETSc CN nonlinear pendulum.

The example itself creates figures; this verifier only opens its
actual files, checks every original accepted GetFEM coefficient,
both MPI rank trajectories, 12-second forced/unforced solutions,
endpoint CN differential and power residuals, physical Hamiltonian
chain-rule accounting, RK4 reference and dt=.01/.005 refinement.
"""
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
from PIL import Image


SCENARIOS=(
    ("preview_free","free",0.,1,.08,.01,True),
    ("preview_forced","forced",.1,1,.08,.01,True),
    ("free_mpi1","free",0.,1,12.,.01,True),
    ("free_mpi2","free",0.,2,12.,.01,False),
    ("free_refined","free",0.,1,12.,.005,False),
    ("forced_mpi1","forced",.1,1,12.,.01,True),
    ("forced_mpi2","forced",.1,2,12.,.01,False),
    ("forced_refined","forced",.1,1,12.,.005,False),
)


def main():
    root=Path(sys.argv[1]).resolve()
    errors=[]
    reports={}
    traces={}
    figures={}
    to_zip=set()
    for label,variant,amplitude,mpi,tf,dt,showcase in SCENARIOS:
        code=root/f"{label}.exitcode"
        if not code.is_file() or code.read_text().strip()!="0":
            errors.append("Original CN run did not complete: "+label)
            continue
        count=round(tf/dt)+1
        for rank in range(mpi):
            prefix=root/f"native_{label}"/f"rank{rank}"
            meta_path=Path(str(prefix)+".json")
            npz_path=Path(str(prefix)+".npz")
            key=(label,rank)
            try:
                meta=json.loads(meta_path.read_text())
                with np.load(npz_path,allow_pickle=False) as source:
                    data={x:np.asarray(source[x]) for x in source.files}
                reports[key]=meta
                traces[key]=data
                to_zip.update((meta_path,npz_path))
                expected_fields=(
                    "theta","p","e_theta","e_p","H",
                    "torque","input_power","viscous_dissipation",
                    "RK4_reference_theta","RK4_reference_p",
                )
                step_fields=(
                    "CN_ODE_theta_defect","CN_ODE_p_defect",
                    "CN_coenergy_rate","CN_step_dissipation",
                    "CN_step_input_power","CN_coenergy_power_defect",
                    "CN_nonlinear_chain_rule_defect",
                    "CN_chain_accounting_defect",
                )
                valid=(
                    meta.get("schema")==
                        "scrimp-v130-original-nonlinear-pendulum-CN-physics-v1"
                    and meta.get("status")==
                        "PASS_ORIGINAL_NONLINEAR_PENDULUM_CN_NATIVE"
                    and meta.get("time_scheme")==
                        "PETSc_TSCN_endpoint_trapezoidal"
                    and meta.get("initialization")=="exact_constitutive"
                    and meta.get("pseudo_TS_executed") is False
                    and meta.get("original_fully_nonlinear_sine") is True
                    and meta.get("canonical_showcase") is showcase
                    and abs(meta.get("torque_amplitude",99.)-amplitude)<1e-14
                    and meta.get("DAE_DOFs")==4
                    and meta.get("steps")==count-1
                    and meta.get("saved_states")==count
                    and meta.get("gates")
                    and all(meta["gates"].values())
                    and data["z"].shape==(count,4)
                    and data["t"].shape==(count,)
                    and all(data[n].shape==(count,) for n in expected_fields)
                    and all(data[n].shape==(count-1,) for n in step_fields)
                    and data["CN_discrete_work"].shape==(count,)
                    and data["CN_cumulative_H_minus_work"].shape==(count,)
                    and np.all(np.isfinite(data["z"]))
                    and np.max(np.abs(data["t"]-np.arange(count)*dt))<1e-9
                )
                if not valid:
                    errors.append(f"Original GetFEM CN time/physics gates {key}")
                    continue
                th,p,et,ep,H=(data[x] for x in
                              ("theta","p","e_theta","e_p","H"))
                t=data["t"]
                tau=amplitude*np.sin(2*np.pi*.3*t)
                m_ep=.5*(ep[1:]+ep[:-1])
                m_et=.5*(et[1:]+et[:-1])
                m_tau=.5*(tau[1:]+tau[:-1])
                Dth=np.diff(th)/dt
                Dp=np.diff(p)/dt
                coenergy=m_et*Dth+m_ep*Dp
                D=.08*m_ep**2
                P=m_tau*m_ep
                nonlinear_chain=np.diff(H)/dt-coenergy
                # Independent algebra, differing grouping only by roundoff.
                calculated={
                    "CN_ODE_theta_defect":Dth-m_ep,
                    "CN_ODE_p_defect":Dp+m_et+.08*m_ep-m_tau,
                    "CN_coenergy_rate":coenergy,
                    "CN_step_dissipation":D,
                    "CN_step_input_power":P,
                    "CN_coenergy_power_defect":coenergy+D-P,
                    "CN_nonlinear_chain_rule_defect":nonlinear_chain,
                }
                for field,value in calculated.items():
                    if not np.allclose(
                        data[field],value,rtol=0.,atol=1e-11
                    ):
                        errors.append(f"Non-native CN step array {key}/{field}")
                if (max(abs(th[0]-1.2),abs(p[0]),
                        abs(et[0]-9.81*np.sin(1.2)),abs(ep[0]))>1e-12
                        or np.max(abs(et-9.81*np.sin(th)))>1e-8
                        or np.max(abs(ep-p))>1e-8
                        or np.max(abs(H-(.5*p*p+9.81*(1-np.cos(th)))))>1e-8
                        or np.max(abs(Dth-m_ep))>1e-7
                        or np.max(abs(Dp+m_et+.08*m_ep-m_tau))>1e-7
                        or np.max(abs(coenergy+D-P))>1e-7):
                    errors.append(f"Original nonlinear CN physical law {key}")
                accumulated=np.r_[0.,np.cumsum(dt*(P-D))]
                if (np.max(abs(accumulated-data["CN_discrete_work"]))>1e-10
                        or np.max(abs(H-H[0]-accumulated
                                      -data["CN_cumulative_H_minus_work"]))>1e-10
                        or np.max(abs(
                            data["CN_cumulative_H_minus_work"][1:]
                            -np.cumsum(dt*nonlinear_chain)
                        ))>1e-8):
                    errors.append(f"Nonquadratic CN chain ledger {key}")
                if amplitude==0 and np.max(np.diff(H))>1e-7:
                    errors.append(f"Free damped CN energy unexpectedly rose {key}")
            except Exception as exc:
                errors.append(f"Original CN evidence unreadable {key}: {exc}")

        if showcase:
            folder=root/f"work_{label}"/"outputs"/f"nonlinear_pendulum_{variant}"
            manifest=folder/"json"/"manifest.json"
            history=folder/"npz"/"native_lumped_time_history.npz"
            try:
                info=json.loads(manifest.read_text())
                to_zip.update((manifest,history))
                if (info.get("status")!=
                        "EXAMPLE_NATIVE_GETFEM_NONLINEAR_PENDULUM_VISUALS"
                        or info.get("time_stepper")!="cn"
                        or info.get("original_saved_GetFEM_DAE_states")!=count
                        or info.get("original_DAE_DOFs")!=4
                        or info.get("variant")!=variant
                        or info.get("torque_amplitude")!=amplitude
                        or info.get("GIF_genuine_accepted_frames")!=min(36,count)
                        or info.get("no_interpolated_or_spatial_FE_states") is not True
                        or info.get("motion_is_schematic") is not True
                        or info.get("CN_is_trapezoidal_not_implicit_midpoint") is not True
                        or info.get("nonquadratic_H_CN_chain_rule_defect_not_assumed_zero") is not True):
                    errors.append(f"Original canonical CN figure manifest {label}")
                if reports.get((label,0),{}).get("canonical_visuals")!=info:
                    errors.append(f"Original canonical example manifest inconsistency {label}")
                with np.load(history,allow_pickle=False) as f:
                    saved={name:np.asarray(f[name]) for name in f.files}
                if (label,0) in traces:
                    original=traces[label,0]
                    direct=(
                        ("t","t"),("theta","theta"),("p","p"),
                        ("e_theta","e_theta"),("e_p","e_p"),("H","H"),
                    )
                    for x,y in direct:
                        if not np.array_equal(saved[x],original[y]):
                            errors.append(f"Canonical CN plot not true GetFEM {label}/{x}")
                    derived=(
                        ("torque","torque"),("input_power","input_power"),
                        ("viscous_dissipation","viscous_dissipation"),
                        ("CN_coenergy_rate","CN_coenergy_rate"),
                        ("CN_power_defect","CN_coenergy_power_defect"),
                        ("CN_step_input_power","CN_step_input_power"),
                        ("CN_step_dissipation","CN_step_dissipation"),
                        ("CN_chain_rule_defect","CN_nonlinear_chain_rule_defect"),
                        ("CN_discrete_work","CN_discrete_work"),
                        ("CN_cumulative_ledger_difference",
                         "CN_cumulative_H_minus_work"),
                    )
                    for x,y in derived:
                        if not np.allclose(
                            saved[x],original[y],rtol=0.,atol=1e-10
                        ):
                            errors.append(f"Canonical CN figure derived values {label}/{x}")
                names=info.get("figures",[])
                if (len(names)!=7 or len(set(names))!=7
                        or "pendulum_motion.gif" not in names):
                    errors.append(f"Original canonical figure inventory {label}")
                for name in names:
                    file=folder/name.rsplit(".",1)[-1]/name
                    if not file.is_file() or file.stat().st_size<=1200:
                        errors.append(f"Missing original canonical figure {label}/{name}")
                    else:
                        to_zip.add(file)
                        if name.endswith(".png"):
                            with Image.open(file) as im:
                                im.verify()
                gif=folder/"gif"/"pendulum_motion.gif"
                with Image.open(gif) as im:
                    if im.n_frames!=min(36,count) or min(im.size)<200:
                        errors.append(f"Not genuine state-count canonical GIF {label}")
                figures[label]={
                    "native_original_states":count,
                    "canonical_path":str(folder.relative_to(root)),
                    "original_PNG_SVG_GIF_files":names,
                    "GIF_accepted_GetFEM_frames":min(36,count),
                }
            except Exception as exc:
                errors.append(f"Original user-example CN figures invalid {label}: {exc}")

    mpi_parity={}
    refinement={}
    for variant in ("free","forced"):
        coarse=traces.get((f"{variant}_mpi1",0))
        finer=traces.get((f"{variant}_refined",0))
        if coarse is None or finer is None:
            errors.append("Missing original CN refinement for "+variant)
        else:
            if (not np.array_equal(coarse["t"],finer["t"][::2])
                    or not np.array_equal(coarse["RK4_reference_theta"],
                                          finer["RK4_reference_theta"][::2])):
                # Each independently integrated high-accuracy RK4 path
                # may disagree at roundoff; time points MUST coincide.
                if not np.allclose(
                    coarse["t"],finer["t"][::2],atol=1e-10,rtol=0.
                ):
                    errors.append(f"Refinement time grid not nested {variant}")
            X=np.column_stack((coarse["theta"],coarse["p"]))
            Y=np.column_stack((finer["theta"],finer["p"]))
            Ref=np.column_stack((
                coarse["RK4_reference_theta"],coarse["RK4_reference_p"]
            ))
            RefFine=np.column_stack((
                finer["RK4_reference_theta"],finer["RK4_reference_p"]
            ))
            error_coarse=float(np.max(np.abs(X-Ref)))
            error_fine=float(np.max(np.abs(Y-RefFine)))
            ratio=error_coarse/error_fine if error_fine>0. else 0.
            refinement[variant]={
                "dt":.01,"refined_dt":.005,
                "coarse_max_abs_theta_or_p_error":error_coarse,
                "refined_max_abs_theta_or_p_error":error_fine,
                "error_ratio":ratio,
                "observed_order_log2_ratio":float(np.log2(ratio))
                    if ratio>0. else None,
                "full_states_coarse":len(X),
                "full_states_refined":len(Y),
            }
            if (error_coarse>=.02 or error_fine>=.006
                    or not (3.3<=ratio<=4.7)):
                errors.append(f"Original CN 12s O(dt²) refinement failed {variant}: "
                              f"{error_coarse}/{error_fine} ratio={ratio}")

        mpi1=traces.get((f"{variant}_mpi1",0))
        mpi20=traces.get((f"{variant}_mpi2",0))
        mpi21=traces.get((f"{variant}_mpi2",1))
        if not all(x is not None for x in (mpi1,mpi20,mpi21)):
            errors.append("Missing full 12-second CN MPI parity "+variant)
            continue
        for rank,other in enumerate((mpi20,mpi21)):
            diffs={}
            for name in (
                "t","z","theta","p","e_theta","e_p","H",
                "torque","input_power","viscous_dissipation",
                "CN_ODE_theta_defect","CN_ODE_p_defect",
                "CN_coenergy_power_defect","CN_nonlinear_chain_rule_defect",
                "CN_discrete_work","CN_cumulative_H_minus_work",
            ):
                if mpi1[name].shape!=other[name].shape:
                    errors.append(f"Native CN MPI shape mismatch {variant}/{name}")
                    continue
                diff=float(np.max(abs(mpi1[name]-other[name])))
                diffs[name]=diff
                if not np.isfinite(diff) or diff>1e-8:
                    errors.append(f"Native CN MPI parity {variant}/{name}: {diff}")
            mpi_parity[f"{variant}_MPI1_to_MPI2_rank{rank}"]=diffs
        if not np.array_equal(mpi20["z"],mpi21["z"]):
            errors.append(f"Native MPI2 replicas differ for {variant}")

    summary={
        "schema":"scrimp-v130-original-nonlinear-pendulum-CN-12s-visual-refinement-v1",
        "status":(
            "PASS_ORIGINAL_PENDULUM_CN_12S_MPI1_MPI2_REFINEMENT_AND_VISUALS"
            if not errors else "FAIL_OR_UNQUALIFIED_PENDULUM_CN_12S"
        ),
        "errors":errors,
        "canonical_example":"examples/nonlinear_pendulum.py",
        "time_stepper":"PETSc_TSCN_endpoint_trapezoidal",
        "exact_constitutive_start":True,
        "full_sine_not_linearized":True,
        "real_original_DAE_DOFS":4,
        "total_accepted_states_per_full_coarse_run":1201,
        "total_accepted_states_per_full_refined_run":2401,
        "complete_full_12s_CN_MPI1_MPI2_parity":mpi_parity,
        "independent_full_sine_RK4_refinement":refinement,
        "original_canonical_example_figures":figures,
        "original_nonquadratic_chain_rule_defect_not_set_to_zero":True,
        "no_spatial_FE_or_interpolated_movie_states":True,
        "no_BDF2_startup_assumed_or_used":True,
    }
    summaryfile=root/"summary.json"
    summaryfile.write_text(json.dumps(summary,indent=2)+"\n")
    to_zip.add(summaryfile)
    to_zip.add(root/"manifest.json")
    to_zip.update(root.glob("*.host.log"))
    to_zip.update(root.glob("*.exitcode"))
    # A scientifically rejected plot must not erase its genuine
    # GetFEM trajectory. Preserve the original example's failure
    # diagnostic in the ZIP for numerical post-mortem.
    to_zip.update(root.glob(
        "work_*/outputs/nonlinear_pendulum_*/"
        "*/failed_native_constitutive_trace.npz"
    ))
    archive=Path(str(root)+".zip")
    with zipfile.ZipFile(
        archive,"x",compression=zipfile.ZIP_DEFLATED,
        compresslevel=4,allowZip64=True,
    ) as z:
        for file in sorted(to_zip):
            if file.is_file():
                z.write(file,str(file.relative_to(root)))
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:
            raise RuntimeError("CN pendulum full archive CRC failure")
    print(json.dumps({
        **summary,
        "archive":str(archive),
        "evidence":str(root),
        "ZIP_CRC_pass":True,
    },indent=2),flush=True)
    raise SystemExit(bool(errors))


if __name__=="__main__":
    main()
