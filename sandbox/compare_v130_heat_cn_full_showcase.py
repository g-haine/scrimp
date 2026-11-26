"""Host-only audit of the real 2-D heat CN long showcase and MPI parity.

Reads original GetFEM/PETSc output files from the canonical example.
Does not solve PDE, alter a FE state, or manufacture graphic data.
"""
import hashlib
import json
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np


def main():
    root=Path(sys.argv[1]).resolve()
    errors=[]
    ranks={}
    datasets={}
    short=[]
    for label,n,steps in (
        ("preview",1,4),("mpi1",1,500),("mpi2",2,500)
    ):
        code=root/f"{label}.exitcode"
        if not code.is_file() or code.read_text().strip()!="0":
            errors.append(f"Native {label} did not complete: exitcode missing or nonzero")
            continue
        canonical=root/f"work_{label}/outputs/png/Hamiltonian_Heat_2D.png"
        if not canonical.is_file() or canonical.stat().st_size<1000:
            errors.append(f"Canonical heat.py Hamiltonian PNG missing: {label}")
        else:
            short.append(canonical)
        for rank in range(n):
            key=(label,rank)
            base=root/f"native_{label}/rank{rank}"
            meta=Path(str(base)+".json")
            data=Path(str(base)+".npz")
            if not meta.is_file() or not data.is_file():
                errors.append(f"Missing original GetFEM FE evidence: {key}")
                continue
            try:
                report=json.loads(meta.read_text(encoding="utf-8"))
                with np.load(data) as archive:
                    z={k:archive[k] for k in archive.files}
                ranks[key]=report
                datasets[key]=z
                init=report.get("initialization_report") or {}
                phys=(report.get("physics") or {})
                gates=phys.get("gates") or {}
                selected=report.get("init_profile_applied") or {}
                if (report.get("status")!="PASS_CN_HEAT_NATIVE_PHYSICS"
                        or report.get("time_stepper")!="cn"
                        or report.get("initializer")!="algebraic"
                        or report.get("init_profile")!="strict"
                        or report.get("DOFs")!=5830
                        or report.get("steps")!=steps
                        or init.get("status")!="PASS"
                        or init.get("accuracy_target_met") is not True
                        or init.get("differential_initial_exact") is not True
                        or init.get("pseudo_ts_executed") is not False
                        or float(init.get("accuracy_target",1.))!=1e-13
                        or selected.get("init_algebraic_max_passes")!="8"
                        or selected.get("init_algebraic_accuracy")!="1e-13"
                        or selected.get("init_algebraic_ksp_rtol")!="1e-14"
                        or selected.get("init_algebraic_ksp_atol")!="1e-14"
                        or not gates or not all(gates.values())
                        or phys.get("steps")!=steps):
                    errors.append(f"Authentic original CN/DAE/energy gate: {key}")
                t=z["time"]; states=z["states"]; H=z["H"]
                if (t.shape!=(steps+1,)
                        or states.shape!=(steps+1,5830)
                        or H.shape!=(steps+1,)
                        or not np.all(np.isfinite(states))
                        or not np.all(np.isfinite(H))
                        or not np.allclose(
                            t,np.arange(steps+1)*0.01,
                            atol=1e-9,rtol=0
                        )):
                    errors.append(f"FE shape, genuine clock or finite values: {key}")
                d=z["Fourier_dissipation"]
                defect=z["CN_step_defect"]
                if (d.shape!=t.shape or defect.shape!=t.shape
                        or not np.all(np.isfinite(d))
                        or not np.all(np.isfinite(defect))
                        or np.min(d)<-1e-10
                        or np.max(np.abs(defect))>1e-8):
                    errors.append(f"Original Fourier or exact CN step power: {key}")
                if not np.all(np.diff(H) < 1e-8):
                    errors.append(f"Unexpected heat Hamiltonian increase: {key}")
            except Exception as exc:
                errors.append(f"Unreadable real FE evidence {key}: {type(exc).__name__}: {exc}")

    parity={}
    if all(k in datasets for k in (("mpi1",0),("mpi2",0),("mpi2",1))):
        reference=datasets["mpi1",0]
        for key in (("mpi2",0),("mpi2",1)):
            other=datasets[key]
            for name in ("time","states","H","Fourier_dissipation",
                         "power_0","power_1","power_2","power_3"):
                if reference[name].shape!=other[name].shape:
                    errors.append(f"Wrong MPI FE or physical power shape: {name} {key}")
            try:
                diff=float(np.max(np.abs(
                    reference["states"]-other["states"]
                )))
                h=float(np.max(np.abs(reference["H"]-other["H"])))
                d=float(np.max(np.abs(
                    reference["Fourier_dissipation"]
                    -other["Fourier_dissipation"]
                )))
                p=float(max(np.max(np.abs(
                    reference[f"power_{k}"]-other[f"power_{k}"]
                )) for k in range(4)))
                parity[f"MPI1_vs_{key[0]}_rank{key[1]}"]={
                    "max_complete_FE_coefficient_abs_diff":diff,
                    "max_H_abs_diff":h,
                    "max_Fourier_dissipation_abs_diff":d,
                    "max_boundary_control_power_abs_diff":p,
                }
                if diff>1e-8 or h>1e-9 or d>1e-9 or p>1e-9:
                    errors.append(f"Full time-history MPI parity failed: {key}")
            except Exception as exc:
                errors.append(f"Cannot measure full MPI parity {key}: {exc}")
        if not np.array_equal(
            datasets["mpi2",0]["states"],datasets["mpi2",1]["states"]
        ):
            errors.append("MPI2 rank replicas are not byte-identical")

    # This actually exercises the canonical heat.py user visualization
    # path, not a one-off sandbox renderer.
    visual_summary={}
    sampled_vtu=[]
    for label,steps in (("preview",4),("mpi1",500)):
        base=root/f"work_{label}/outputs"
        figures=base/"showcase/heat/figures"
        manifest=base/"showcase/heat/manifest.json"
        vtu_collection=base/"pv/T/T.pvd"
        original=ranks.get((label,0),{})
        if not manifest.is_file() or not vtu_collection.is_file():
            errors.append(f"Missing canonical original heat showcase or PVD {label}")
            continue
        try:
            metadata=json.loads(manifest.read_text(encoding="utf-8"))
            if (metadata.get("visual_status")!="GENERATED_FROM_ACTUAL_GETFEM_NATIVE_FE"
                    or metadata.get("VTK_frames")!=steps+1
                    or metadata.get("GIF_original_FE_frames")!=min(steps+1,36)
                    or metadata.get("no_spatial_or_time_interpolation") is not True):
                errors.append(f"Real-FE showcase manifest assertions failed: {label}")
            if (original.get("visuals") or {}).get("visual_status")!=(
                    "GENERATED_FROM_ACTUAL_GETFEM_NATIVE_FE"):
                errors.append(f"Original heat.py rank0 did not publish visuals {label}")
            expected=(
                "heat_temperature_snapshots.png",
                "heat_temperature_snapshots.svg",
                "heat_energy_and_ports.png",
                "heat_energy_and_ports.svg",
                "heat_temperature_evolution.gif",
            )
            for name in expected:
                p=figures/name
                if not p.is_file() or p.stat().st_size<1000:
                    errors.append(f"Missing/empty original heat visual: {label}/{name}")
                else:
                    short.append(p)
            # Check actual GIF frames, not only the name or metadata.
            from PIL import Image
            with Image.open(figures/"heat_temperature_evolution.gif") as gif:
                if gif.n_frames!=min(steps+1,36) or min(gif.size)<200:
                    errors.append(f"Invalid native heat animation frames: {label}")
            short.extend((manifest,vtu_collection))
            collections=ET.parse(vtu_collection).getroot().findall(".//DataSet")
            if len(collections)!=steps+1:
                errors.append(f"Original heat VTK clock length mismatch: {label}")
                continue
            time=datasets[label,0]["time"]
            inventory=[]
            for i,node in enumerate(collections):
                filename=Path(node.attrib.get("file",""))
                if (filename.name!=str(filename)
                        or filename.suffix!=".vtu"):
                    errors.append(f"Unsafe canonical VTK filename at {label}/{i}")
                    break
                if abs(float(node.attrib["timestep"])-time[i])>1e-9:
                    errors.append(f"Original heat VTK timestamp mismatch {label}/{i}")
                    break
                f=vtu_collection.parent/filename
                if not f.is_file() or f.stat().st_size<100:
                    errors.append(f"Missing native VTK at {label}/{i}")
                    break
                # Hash ALL 501 original VTK files for archival evidence;
                # retain first/third/middle/last raw VTUs in compact ZIP.
                digest=hashlib.sha256(f.read_bytes()).hexdigest()
                inventory.append({
                    "index":i,"time":float(time[i]),"filename":filename.name,
                    "bytes":f.stat().st_size,"SHA256":digest,
                })
                if i in (0,steps//3,2*steps//3,steps):
                    sampled_vtu.append(f)
            if len(inventory)!=steps+1:
                errors.append(f"Incomplete canonical original VTK inventory {label}")
            inv=base/"showcase/heat/native_vtk_inventory.json"
            inv.write_text(json.dumps({
                "schema":"scrimp-heat-native-vtk-original-full-inventory-v1",
                "complete_native_frame_count":len(inventory),
                "all_original_VTK_on_host":True,
                "ZIP_contains_selected_raw_VTK_only":True,
                "frames":inventory,
            },indent=2)+"\n")
            short.append(inv)
            visual_summary[label]={
                "full_VTK_frame_count":len(inventory),
                "original_GetFEM_snapshots_in_GIF":min(steps+1,36),
                "figures":[str((figures/x).relative_to(root)) for x in expected],
                "native_VTK_manifest":str(inv.relative_to(root)),
            }
        except Exception as exc:
            errors.append(f"Original heat.py native showcase failed {label}: "
                          f"{type(exc).__name__}: {exc}")

    evidence={
        "schema":"scrimp-v130-heat-CN-strict-full-physical-and-visual-v1",
        "status":("PASS_NATIVE_HEAT_CN_500_STEPS_MPI1_MPI2_AND_VISUALS"
                  if not errors else "FAIL_OR_UNQUALIFIED_HEAT_CN_LONG"),
        "errors":errors,
        "original_example":"examples/heat.py",
        "init_profile":"strict_by_default_for_cn",
        "physical_scheme":"fixed_CN_0.01",
        "full_real_FE_states":501,
        "original_mesh_size":0.1,
        "original_dof_count":5830,
        "full_FE_MPI_parity":parity,
        "native_qualifications":{
            f"{label}_rank{rank}":{
                "steps":val.get("steps"),
                "initialization_residual":(
                    val.get("initialization_report") or {}
                ).get("native_final_residual_l2"),
                "max_abs_CN_step_energy_defect":(
                    val.get("physics") or {}
                ).get("max_abs_CN_step_energy_defect"),
                "max_abs_CN_cumulative_energy_defect":(
                    val.get("physics") or {}
                ).get("max_abs_CN_cumulative_energy_defect"),
                "H_initial":(val.get("physics") or {}).get("H_initial"),
                "H_final":(val.get("physics") or {}).get("H_final"),
            } for (label,rank),val in ranks.items()
        },
        "visualization":visual_summary,
        "complete_original_VTK_directory_on_user_host":True,
        "ZIP_contains_only_hashed_selected_VTK_snapshots":True,
    }
    summary=root/"summary.json"
    summary.write_text(json.dumps(evidence,indent=2)+"\n")
    archive_files=set(short)
    archive_files.add(summary)
    for f in root.glob("*.exitcode"):
        archive_files.add(f)
    for f in root.glob("*.host.log"):
        archive_files.add(f)
    for f in root.glob("*.json"):
        archive_files.add(f)
    for folder in ("native_preview","native_mpi1","native_mpi2"):
        archive_files.update(
            f for f in (root/folder).rglob("*") if f.is_file()
        )
    inclusion=root/"archive_inclusion_manifest.json"
    inclusion.write_text(json.dumps({
        "schema":"scrimp-heat-compact-evidence-manifest-v1",
        "full_original_VTK_retained_on_user_host":True,
        "complete_original_VTK_hashed_in_showcase_inventory":True,
        "included_files":sorted(str(f.relative_to(root)) for f in archive_files),
    },indent=2)+"\n")
    archive_files.add(inclusion)
    import zipfile
    archive=Path(str(root)+".zip")
    with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,
                         compresslevel=4,allowZip64=True) as z:
        for f in sorted(archive_files):
            if f.is_file():z.write(f,arcname=str(f.relative_to(root)))
    print(json.dumps({
        **evidence,"archive":str(archive),"evidence":str(root)
    },indent=2),flush=True)
    raise SystemExit(bool(errors))


if __name__=="__main__":
    main()
