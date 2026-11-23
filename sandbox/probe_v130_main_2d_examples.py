# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native 2-D main-example diagnostic pilot: original physics + output files.

Run **MPI=1 only**: both existing top-level functions include ParaView
export postprocessing, which is NOT certified for concurrent MPI writers.
This is a real GetFEM/PETSc/Python integration test, not a performance
benchmark or a claim of universal BDF discrete energy conservation.
The heat-wave and SWE models are called without modifying weak forms.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

from examples.heat_wave import heat_wave_eq
from examples.shallow_water import shallow_water_eq
from measure_v130_heat_wave_physics import collect_heat_wave_physics


T_FINAL = 0.004
DT = 0.001
BACKENDS = ("mumps", "superlu_dist")


def real_array(values, name):
    array=np.asarray(values)
    if np.iscomplexobj(array):
        imag=float(np.max(abs(array.imag)))
        if imag>1e-10:
            raise RuntimeError(f"Unexpected complex {name}: imag {imag}")
        array=array.real
    array=np.asarray(array,dtype=float)
    if not np.all(np.isfinite(array)):
        raise RuntimeError("Nonfinite original "+name)
    return array


def check_output_files(example):
    """Original public plot and ParaView outputs; no synthetic replacement."""
    base=Path.cwd()/"outputs"
    image=(
        "Hamiltonian_Heat1_Wave2_2D.png"
        if example=="heat_wave" else
        "Hamiltonian_Inviscid_Shallow_Water_2D.png"
    )
    plot=base/"png"/image
    if not plot.is_file() or plot.stat().st_size<1000:
        raise RuntimeError(f"Original public Hamiltonian PNG missing: {plot}")
    variables=("p","T") if example=="heat_wave" else ("h",)
    result={}
    for var in variables:
        folder=base/"pv"/var
        pvd=folder/f"{var}.pvd"
        vtus=sorted(folder.glob("*.vtu"))
        if not pvd.is_file() or pvd.stat().st_size<100:
            raise RuntimeError(f"Original public ParaView PVD missing: {pvd}")
        if not vtus or any(p.stat().st_size<100 for p in vtus):
            raise RuntimeError(f"Original public ParaView VTU data missing: {folder}")
        result[var]={"pvd":str(pvd.relative_to(base)),
                     "vtu_count":len(vtus),
                     "first_vtu_bytes":vtus[0].stat().st_size}
    return image, result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--example",choices=("heat_wave","shallow_water"),
                   required=True)
    p.add_argument("--backend",choices=BACKENDS,required=True)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    comm=PETSc.COMM_WORLD.tompi4py()
    if comm.size!=1:
        raise ValueError(
            "Original 2D example ParaView exports are not yet MPI-safe; "
            "pilot requires MPI=1"
        )
    if not PETSc.Sys.hasExternalPackage(args.backend):
        raise RuntimeError("PETSc LU backend absent: "+args.backend)

    # Invoke actual shipped top-level example including original
    # plot_Hamiltonian and export_to_pv postprocessing.
    if args.example=="heat_wave":
        system=heat_wave_eq(
            heat_region=1,wave_region=2,final_time=T_FINAL,
            time_step=DT,save_step=DT,
            factor_solver=args.backend,
        )
        names=("T","p","q")
    else:
        system=shallow_water_eq(
            final_time=T_FINAL,time_step=DT,save_step=DT,
            factor_solver=args.backend,
        )
        names=("h","p")

    chosen=system.time_scheme.getString(
        "pc_factor_mat_solver_type",""
    ).lower()
    if chosen!=args.backend:
        raise RuntimeError(f"Original example selected {chosen} != {args.backend}")
    times=real_array(system.solution["t"],"physical times")
    states=real_array(system.solution["z"],"complete FE state")
    H=real_array(system.get_Hamiltonian(),"Hamiltonian")
    # The original BDF2 example is ADAPTIVE. dt is initial/proposed
    # step size, not a promise of exactly four accepted steps.
    # Likewise the normal SCRIMP monitor saves approximately every
    # dt_save but does not interpolate to an exact uniform clock.
    # The previous equality-to-four gate incorrectly rejected a
    # physically completed 39-accepted, 4-rejected PETSc integration.
    if (times.ndim!=1 or times.size<2
            or not np.all(np.isfinite(times))
            or abs(float(times[0]))>1e-12
            or not np.all(np.diff(times)>0.)
            or states.ndim!=2 or states.shape[0]!=times.size
            or H.shape!=(times.size,)
            or not np.all(np.isfinite(states))
            or system.last_ts_converged_reason<=0
            or system.last_ts_rejections<0
            or system.last_ts_snes_failures!=0
            or system.last_ts_steps < times.size-1
            or abs(float(times[-1])-T_FINAL)>1e-11
            or abs(system.last_ts_time-T_FINAL)>1e-11):
        raise RuntimeError("Original adaptive 2D PETSc trajectory/TS failure")
    if np.min(H)<=0.:
        raise RuntimeError("Original physical Hamiltonian is not positive")
    original_vars={}
    for name in names:
        x=real_array(system.get_solution(name),"original field "+name)
        if x.ndim!=2 or x.shape[0]!=times.size:
            raise RuntimeError("Incorrect original full field "+name)
        original_vars[name]=x
    if float(np.max(np.linalg.norm(states-states[0],axis=1)))<1e-10:
        raise RuntimeError("Original PDE dynamics appear completely inactive")
    sample_mass=None
    if args.example=="shallow_water":
        sample_mass=real_array(
            system.get_quantity("h",region=1),
            "original SWE integrated water height"
        )
        if sample_mass.shape!=(times.size,) or np.min(sample_mass)<=0.:
            raise RuntimeError("Invalid original SWE water mass")
        # Preserve the observed mass defect for model review. Do not
        # confuse a BDF2 energy/time-discretization error with LU parity.
    figure,paraview=check_output_files(args.example)

    # Keep the heat-wave original-physics audit SEPARATE from
    # comparing the factorization backends. A solver PASS must
    # not be mistaken for an interface/power-balance PASS.
    physics_data={}
    physics_report=None
    if args.example=="heat_wave":
        physics_data,physics_report=collect_heat_wave_physics(
            system,times,H
        )
    args.out.mkdir(parents=True,exist_ok=True)
    trace=args.out/"original_2d_example_trace.npz"
    report=args.out/"original_2d_example_report.json"
    if trace.exists() or report.exists():
        raise FileExistsError("Refuse to overwrite existing native evidence")
    np.savez_compressed(
        trace,times=times,states=states,H=H,
        **original_vars,
        **({"integral_h":sample_mass} if sample_mass is not None else {}),
        **({"phys_"+name:value
           for name,value in physics_data.items()} if physics_data else {})
    )
    if physics_report is not None:
        (args.out/"heat_wave_physics_diagnostic.json").write_text(
            json.dumps(physics_report,indent=2)+"\n",
            encoding="utf-8"
        )
    result={
        "schema":"scrimp-v130-original-main-2d-example-pilot-v1",
        "status":(
            "PASS_NATIVE_INTEGRATION_HEAT_WAVE_PHYSICAL_AUDIT_PENDING"
            if args.example=="heat_wave"
            else "PASS_NATIVE_ORIGINAL_2D_EXAMPLE_SINGLE_MPI"
        ),
        "example":args.example,"factor_solver":chosen,
        "original_FEM_and_weak_forms_unchanged":True,
        "PETSc_reason":int(system.last_ts_converged_reason),
        "rejected_steps":int(system.last_ts_rejections),
        "SNES_failures":int(system.last_ts_snes_failures),
        "steps":int(system.last_ts_steps),
        "requested_initial_dt":DT,
        "actual_saved_time_count":int(times.size),
        "accepted_min_dt":float(np.min(np.diff(times))),
        "accepted_max_saved_dt":float(np.max(np.diff(times))),
        "saved_time_grid_is_uniform":bool(
            np.allclose(np.diff(times),np.diff(times)[0],
                        atol=1e-12,rtol=0.)
        ),
        "adaptive_integrator_original_behavior_preserved":True,
        "MPI":comm.size,"full_state_DOFs":int(states.shape[1]),
        "t_initial":float(times[0]),"t_final":float(times[-1]),
        "dt":DT,"Hamiltonian_initial":float(H[0]),
        "Hamiltonian_final":float(H[-1]),
        "Hamiltonian_delta":float(H[-1]-H[0]),
        "original_water_height_integral_initial":(
            None if sample_mass is None else float(sample_mass[0])
        ),
        "original_water_height_integral_final":(
            None if sample_mass is None else float(sample_mass[-1])
        ),
        "original_outputs_plot":figure,
        "original_outputs_paraview":paraview,
        "BDF_energy_balance_identity_not_claimed":True,
        "physical_model_approval_not_inferred_from_solver_success":True,
        "parallel_export_safety_not_qualified":True,
        "heat_wave_physical_audit":physics_report,
        "heat_wave_physics_certified":False,
    }
    report.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(
        "PASS ORIGINAL 2D EXAMPLE " + args.example
        + f" LU={chosen} MPI=1 DOFs={states.shape[1]}"
        + f" accepted_steps={system.last_ts_steps}"
        + f" rejected={system.last_ts_rejections}"
        + f" saved={times.size}"
        + f" dH={H[-1]-H[0]:.8g}"
        + (" integral_h_change="
           + f"{sample_mass[-1]-sample_mass[0]:.8g}"
           if sample_mass is not None else ""),
        flush=True,
    )


if __name__=="__main__":
    main()
