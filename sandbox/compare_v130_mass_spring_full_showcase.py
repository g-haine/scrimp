# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Host-only, fail-closed verification of original canonical oscillator.

Authentic outputs are produced by examples/mass_spring_damper.py,
not by this checker. Test all 801 original accepted DAE states,
both MPI2 replicas, true BDF2 identities and source-generated figures.
"""
import json
import sys
import zipfile
from pathlib import Path

import numpy as np


def main():
    root=Path(sys.argv[1]).resolve()
    errors=[]
    ranks={}
    traces={}
    archive_files=set()
    figures={}
    for label,size,count,showcase in (
        ("preview",1,9,True),
        ("mpi1",1,801,True),
        ("mpi2",2,801,False),
    ):
        exitfile=root/f"{label}.exitcode"
        if not exitfile.is_file() or exitfile.read_text().strip()!="0":
            errors.append(f"Original native run not completed: {label}")
            continue
        for rank in range(size):
            key=(label,rank)
            path=root/f"native_{label}/rank{rank}"
            j=Path(str(path)+".json")
            z=Path(str(path)+".npz")
            if not j.is_file() or not z.is_file():
                errors.append(f"Missing original GetFEM full DAE evidence {key}")
                continue
            archive_files.update((j,z))
            try:
                meta=json.loads(j.read_text())
                with np.load(z,allow_pickle=False) as npz:
                    trace={k:np.asarray(npz[k]) for k in npz.files}
                ranks[key]=meta
                traces[key]=trace
                t=trace["t"]
                all_physics=meta.get("gates") or {}
                target=(count-1)*0.01
                if (meta.get("status")!="PASS_ORIGINAL_LUMPED_ODE_NATIVE_PILOT"
                        or meta.get("initializer")!="exact_constitutive"
                        or meta.get("no_pseudo_time_initialization") is not True
                        or meta.get("time_stepper")!="original_BDF2"
                        or meta.get("canonical_showcase") is not showcase
                        or meta.get("steps")!=count-1
                        or meta.get("saved_states")!=count
                        or meta.get("full_DAE_DOFs")!=4
                        or abs(meta.get("t_final",0.)-target)>1e-12
                        or not all_physics or not all(all_physics.values())):
                    errors.append(f"Native original BDF2/ODE gate failed: {key}")
                if (t.shape!=(count,)
                        or trace["z"].shape!=(count,4)
                        or any(trace[x].shape!=(count,) for x in (
                            "q","p","e_q","e_p","H",
                            "external_force","external_power",
                            "Fourier_or_viscous_dissipation",
                            "reference_q","reference_p"
                        ))
                        or trace["BDF2_step_power_defect"].shape!=(count-2,)
                        or not np.all(np.isfinite(trace["z"]))
                        or np.max(np.abs(t-np.arange(count)*.01))>1e-10):
                    errors.append(f"Original FE coefficients/time grid invalid {key}")
                elif (np.max(np.abs(trace["BDF2_step_power_defect"]))>1e-8
                      or np.max(np.abs(trace["e_q"]-5*trace["q"]))>1e-8
                      or np.max(np.abs(trace["e_p"]-trace["p"]/1.5))>1e-8
                      or max(abs(trace["q"][0]-.4),abs(trace["p"][0]))>1e-12):
                    errors.append(f"Original physics identity or exact initial state {key}")
            except Exception as exc:
                errors.append(f"Cannot independently read original native state {key}: {exc}")

        if showcase:
            folder=root/f"work_{label}/outputs/showcase/mass_spring_damper"
            manifest=folder/"manifest.json"
            time_history=folder/"native_lumped_time_history.npz"
            if not manifest.is_file() or not time_history.is_file():
                errors.append(f"Example-native figure manifest/data missing {label}")
                continue
            archive_files.update((manifest,time_history))
            try:
                metadata=json.loads(manifest.read_text())
                if (metadata.get("visual_status")!=
                        "GENERATED_BY_CANONICAL_EXAMPLE_FROM_GETFEM_LUMPED_STATES"
                        or metadata.get("native_GetFEM_saved_states")!=count
                        or metadata.get("native_original_full_DAE_DOFS")!=4
                        or metadata.get("original_accepted_states_used_for_GIF")
                            !=min(count,36)
                        or metadata.get("no_spatial_FEM_or_time_interpolation")
                            is not True
                        or metadata.get("BDF2_startup_energy_identity_not_claimed")
                            is not True):
                    errors.append(f"Canonical example-native manifest mismatch {label}")
                original=ranks.get((label,0),{})
                if original.get("canonical_visuals")!=metadata:
                    errors.append(f"Native original example did not publish same visual metadata {label}")
                with np.load(time_history,allow_pickle=False) as f:
                    vis={k:np.asarray(f[k]) for k in f.files}
                if (label,0) in traces:
                    original_trace=traces[label,0]
                    for name,signal in (
                        ("time","t"),("q","q"),("p","p"),
                        ("e_q","e_q"),("e_p","e_p"),("H","H"),
                        ("force","external_force"),
                        ("dissipation","Fourier_or_viscous_dissipation"),
                        ("input_power","external_power"),
                        ("BDF2_step_defect","BDF2_step_power_defect"),
                    ):
                        # Raw GetFEM values must match bitwise. Power
                        # terms are independently recomputed from
                        # algebraically equivalent p/m and e_p; rounding
                        # at ~1e-17 is expected and not an invented plot.
                        raw=name in ("time","q","p","e_q","e_p","H","force")
                        ok=(np.array_equal(vis[name],original_trace[signal])
                            if raw else
                            np.allclose(vis[name],original_trace[signal],
                                        atol=1e-12,rtol=0))
                        if not ok:
                            errors.append(f"Not real FE data in canonical figure {label}/{name}")
                for name in metadata.get("figures",[]):
                    file=folder/"figures"/name
                    if not file.is_file() or file.stat().st_size<1200:
                        errors.append(f"Missing original figure or empty plot {label}/{name}")
                    else:
                        archive_files.add(file)
                from PIL import Image
                gif=folder/"figures/mass_spring_motion.gif"
                with Image.open(gif) as img:
                    if img.n_frames!=min(count,36) or min(img.size)<200:
                        errors.append(f"Not genuine accepted-state GIF {label}")
                figures[label]={
                    "native_states":count,
                    "original_GIF_saved_frames":min(count,36),
                    "canonical_folder":str(folder.relative_to(root)),
                    "files":metadata.get("figures"),
                }
            except Exception as exc:
                errors.append(f"Native example visualization inspection failed {label}: {exc}")

    full_parity={}
    if all(k in traces for k in (("mpi1",0),("mpi2",0),("mpi2",1))):
        reference=traces["mpi1",0]
        for key in (("mpi2",0),("mpi2",1)):
            other=traces[key]
            error_by_signal={}
            for name in (
                "t","z","q","p","e_q","e_p","H",
                "external_force","external_power",
                "Fourier_or_viscous_dissipation",
                "BDF2_step_power_defect",
            ):
                if reference[name].shape!=other[name].shape:
                    errors.append(f"Original MPI state shape mismatch {key}/{name}")
                    continue
                diff=float(np.max(np.abs(reference[name]-other[name])))
                error_by_signal[name]=diff
                if not np.isfinite(diff) or diff>1e-8:
                    errors.append(f"Original 801-state MPI parity failed {key}/{name}: {diff}")
            full_parity[f"mpi1_to_{key[0]}_rank{key[1]}"]=error_by_signal
        if not np.array_equal(traces["mpi2",0]["z"],traces["mpi2",1]["z"]):
            errors.append("MPI2 original DAE histories not byte-identical")

    summary={
        "schema":"scrimp-v130-original-mass-spring-8s-mpi-physical-visual-v1",
        "status":("PASS_NATIVE_ORIGINAL_MASS_SPRING_800STEP_MPI1_MPI2_AND_VISUALS"
                  if not errors else "FAIL_OR_UNQUALIFIED_MASS_SPRING_8S"),
        "errors":errors,
        "full_original_time_final":8.,
        "full_original_BDF2_time_step":.01,
        "original_lumped_DAE_DOFs":4,
        "full_accepted_state_count":801,
        "initialization":"exact_constitutive_no_pseudo",
        "initial_BDF2_PETSc_halfstep_identity_not_claimed":True,
        "full_MUMPS_MPI1_MPI2_parity":full_parity,
        "original_example_created_visuals":figures,
        "rank_physics":{
            f"{key[0]}_rank{key[1]}":{
                "native_q0":v.get("native_q0"),
                "native_p0":v.get("native_p0"),
                "max_exact_quadratic_BDF2_step_defect":
                    (v.get("residuals") or {}).get(
                        "exact_quadratic_BDF2_step_energy_k_ge_2"
                    ),
                "RK4_reference_max_q_or_p_error":
                    (v.get("residuals") or {}).get(
                        "RK4_reference_max_q_or_p_error"
                    ),
                "H_initial":v.get("H_initial"),
                "H_final":v.get("H_final"),
            } for key,v in ranks.items()
        },
        "canonical_visuals_are_generated_by_example_not_test":True,
        "no_unqualified_grid_refinement_claim":True,
    }
    summary_path=root/"summary.json"
    summary_path.write_text(json.dumps(summary,indent=2)+"\n")
    archive_files.update((summary_path,root/"manifest.json"))
    archive_files.update(root.glob("*.host.log"))
    archive_files.update(root.glob("*.exitcode"))
    archive=Path(str(root)+".zip")
    with zipfile.ZipFile(
        archive,"x",compression=zipfile.ZIP_DEFLATED,
        compresslevel=4,allowZip64=True
    ) as z:
        for file in sorted(archive_files):
            if file.is_file():
                z.write(file,arcname=str(file.relative_to(root)))
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:
            raise RuntimeError("Scientific original mass-spring ZIP CRC failure")
    print(json.dumps({
        **summary,
        "evidence":str(root),
        "archive":str(archive),
        "ZIP_CRC_pass":True,
    },indent=2),flush=True)
    raise SystemExit(bool(errors))


if __name__=="__main__":
    main()
