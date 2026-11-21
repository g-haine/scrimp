# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.

"""Independent host audit of #28/#6 FULL C1 physical spinup PETSc TSCN.

Check all original full-vector native trajectory samples; real TS metadata,
individual endpoint TSCN and algebraic constraints, h-weighted wall work,
viscous volume dissipation, MPI1/2 coefficient parity, and h/p endpoint
differences for N=10/20/40 with SAME final physical time at end of ramp.

Report temporal differences/order estimates, but DO NOT assert second
order from a single smooth spinup or claim that PETSc endpoint
constant-mass assumptions apply to nonlinear M(h).
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np

OCI = "1d6d0c90115213298047b28237a7f38087563a35"
MESHES = (.2, .1)
STEP_COUNTS = (10, 20, 40)
TF = .02
DOFS = {.2: 3300, .1: 10272}
DYNAMIC_DOFS = {.2: 1452, .1: 4758}


def check(cond, why):
    if not cond:
        raise ValueError(why)


def finite(v, name):
    check(type(v) in (float, int) and math.isfinite(v),
          f"Missing or nonfinite physical {name}: {v!r}")
    return float(v)


def norm(values):
    return float(np.linalg.norm(values))


def validate(data, mesh, ranks, count):
    dt = TF/count
    check(data.get("schema") == "scrimp-v130-issue28-full-spinup-v1"
          and data.get("status") == "NATIVE_TSCN_FULL_C1_RAMP_DIAGNOSTIC",
          "Invalid native physical refinement schema/status")
    check(data.get("nominal_mesh") == mesh
          and data.get("mpi_ranks") == ranks
          and data.get("geometry_order") == 2
          and data.get("height_order") == 2
          and data.get("target_rotation") == .2
          and data.get("initial_height") == 5.
          and data.get("initial_momentum_identically_zero") is True
          and data.get("original_physical_facet_tangent_wall") is True
          and data.get("original_viscous_and_gyro_weak_bricks") is True
          and data.get("original_incompatible_rotating_initial_q0_untouched") is True,
          "Original physical P2 SWE or distinct compatible spinup altered")
    full_control = finite(
        data.get("native_full_target_wall_control_U_l2"), "full rotation wall control"
    )
    plateau = finite(
        data.get("native_post_ramp_plateau_wall_control_U_l2"), "ramp plateau"
    )
    check(data.get("end_of_C1_ramp_reached") is True
          and data.get("original_ramp_start_of_decreasing_envelope") == .03
          and full_control > 1e-12 and plateau > 1e-12
          and math.isclose(full_control,plateau,rel_tol=1e-10),
          "Original physical wall never reached full C1 target rotation")
    check(data.get("physical_scheme") == "petsc_tscn_endpoint_theta_half"
          and data.get("nonlinear_mass_PETSc_endpoint_assumption_violated") is True
          and data.get("nonlinear_mass_time_accuracy_qualified") is False
          and data.get("time_order_qualified") is False
          and data.get("physical_energy_balance_convergence_qualified") is False
          and data.get("multistep_rotating_trajectory_qualified") is False
          and data.get("general_high_index_initializer_qualified") is False,
          "Unjustified TSCN nonlinear-mass/physical-integration qualification")
    init = data.get("initial_analytic_publication")
    check(isinstance(init, dict)
          and init.get("strategy") == "sandbox_exact_stationary_spinup_only"
          and init.get("status") ==
          "SANDBOX_ONLY_PUBLISHED_NOT_GENERAL_INITIALIZER"
          and init.get("differential_q0_preserved_bitwise") is True
          and init.get("getfem_globalvec_published_bitwise") is True
          and init.get("no_pseudo_or_BDF_initialization") is True,
          "Unsafe or unqualified t0 physical-state publication")
    check(finite(init.get("initial_native_IFunction", {}).get("full_l2"),
                 "initial exact DAE full residual") < 1e-10,
          "Initial pH-DAE constraints were not exactly compatible")
    intervals = data.get("dynamic_variable_intervals")
    check(isinstance(intervals, dict) and set(intervals)=={"h","p"},
          "Missing native original GetFEM dynamic intervals")
    selected = []
    for name in ("h","p"):
        iv=intervals[name]
        check(isinstance(iv,list) and len(iv)==2
              and all(type(x) is int for x in iv)
              and 0 <= iv[0] < DOFS[mesh]
              and 0 < iv[1] <= DOFS[mesh]-iv[0],
              f"Invalid native {name} GetFEM interval")
        selected.extend(range(iv[0],iv[0]+iv[1]))
    check(len(set(selected))==DYNAMIC_DOFS[mesh]
          and data.get("actual_total_dofs")==DOFS[mesh]
          and data.get("differential_dofs")==DYNAMIC_DOFS[mesh],
          "Native GetFEM dynamic h/p variable layout changed")
    check(data.get("requested_steps") == count
          and data.get("ts_accepted_steps") == count
          and data.get("ts_rejected_steps") == 0
          and data.get("ts_snes_failures") == 0
          and type(data.get("ts_reason")) is int
          and data.get("ts_reason") > 0
          and math.isclose(finite(data.get("fixed_step"), "dt"), dt,
                           rel_tol=0., abs_tol=1e-14)
          and data.get("final_time_fixed") == TF
          and abs(finite(data.get("actual_TS_time"), "PETSc terminal time")
                  -TF) < 1e-12,
          "Physical TSCN did not reach fixed t_final with requested steps")
    stages = data.get("all_full_native_IFunction_TSCN_step_residuals")
    fa = data.get("all_endpoint_original_algebraic_residuals")
    midpoint = data.get("all_implicit_midpoint_OUT_OF_SCHEME_residuals")
    check(isinstance(stages, list) and len(stages) == count
          and isinstance(fa, list) and len(fa) == count
          and isinstance(midpoint, list) and len(midpoint) == count,
          "Missing per-step original GetFEM TSCN/DAE diagnostics")
    for k, (x, y, z) in enumerate(zip(stages, fa, midpoint)):
        check(isinstance(x, dict), "Non-dictionary TSCN residual")
        for key in ("full_l2", "dynamic_l2", "algebraic_l2"):
            check(0 <= finite(x.get(key), f"TSCN step {k} {key}") < 1e-7,
                  f"Original GetFEM TSCN equation failed at step {k+1}")
        check(0 <= finite(y, f"endpoint algebraic at step {k}") < 1e-7,
              f"Original FE algebraic constraints failed at step {k+1}")
        check(0 <= finite(z, f"out-of-scheme midpoint {k}"),
              "Invalid informational midpoint diagnostic")

    require_fields = (
        "maximum_native_TSCN_equation_residual_l2",
        "maximum_endpoint_native_algebraic_residual_l2",
        "final_dynamic_state_change_l2",
        "minimum_height_coefficient_over_trajectory",
        "initial_hamiltonian", "final_hamiltonian",
        "final_hamiltonian_change", "integrated_physical_wall_work",
        "integrated_viscous_dissipation",
        "integrated_viscous_grad", "integrated_viscous_div",
        "integrated_normal_port_work", "integrated_tangent_port_work",
        "final_physical_power_balance_defect",
        "maximum_physical_power_balance_defect_l2",
        "final_defect_over_actual_energy_exchange",
        "final_volume_balance_defect", "maximum_volume_balance_defect",
    )
    v = {k: finite(data.get(k), k) for k in require_fields}
    check(math.isclose(v["maximum_native_TSCN_equation_residual_l2"],
                       max(x["full_l2"] for x in stages),
                       rel_tol=1e-7, abs_tol=1e-15)
          and math.isclose(v["maximum_endpoint_native_algebraic_residual_l2"],
                           max(fa), rel_tol=1e-7, abs_tol=1e-15),
          "Reported maximum original TSCN/DAE residual differs from actual steps")
    check(v["minimum_height_coefficient_over_trajectory"] > 0
          and v["final_dynamic_state_change_l2"] > 1e-12
          and v["integrated_viscous_dissipation"] >= -1e-12
          and data.get("physical_power_audit_used_h_weighted_boundary_terms") is True,
          "False physical-state evolution/height/dissipation claim")
    deltaH=v["final_hamiltonian_change"]
    W=v["integrated_physical_wall_work"]
    D=v["integrated_viscous_dissipation"]
    defect=deltaH+D-W
    scale=max(abs(deltaH),abs(W),abs(D),1e-14)
    check(math.isclose(defect,v["final_physical_power_balance_defect"],
                       rel_tol=1e-7,abs_tol=2e-12)
          and math.isclose(defect/scale,v["final_defect_over_actual_energy_exchange"],
                           rel_tol=1e-7,abs_tol=1e-6)
          and math.isclose(v["integrated_normal_port_work"]
                           +v["integrated_tangent_port_work"],W,
                           rel_tol=1e-8,abs_tol=1e-12)
          and math.isclose(v["integrated_viscous_grad"]
                           +v["integrated_viscous_div"],D,
                           rel_tol=1e-8,abs_tol=1e-12),
          "Physical weighted work/viscosity/Hamiltonian bookkeeping incorrect")
    return {"Fa_max": v["maximum_endpoint_native_algebraic_residual_l2"],
            "TSCN_max": v["maximum_native_TSCN_equation_residual_l2"],
            "balance_defect": defect,
            "balance_over_exchange": defect/scale,
            "volume_defect": v["final_volume_balance_defect"],
            "delta_H":deltaH, "W":W, "D":D}


def load_native_case(root, mesh, ranks, count):
    path = root/f"result_h{str(mesh).replace('.', 'p')}_n{count}_mpi{ranks}"
    report_path=path/"full_spinup.json"
    npz_path=path/"full_spinup.npz"
    check(report_path.is_file() and npz_path.is_file(),
          f"Missing native TSCN refinement evidence for {mesh=} {count=} {ranks=}")
    report=json.loads(report_path.read_text(encoding="utf-8"))
    metrics=validate(report, mesh, ranks, count)
    with np.load(npz_path, allow_pickle=False) as data:
        required = {"times","states","H","kinetic","potential",
                    "p_normal","p_viscous","d_grad","d_div",
                    "volume","volume_balance","supplied","dissipated",
                    "balance","relative_balance_exchange"}
        check(required <= set(data.files), "Missing original saved native physical traces")
        times=np.asarray(data["times"],dtype=float)
        states=np.asarray(data["states"],dtype=float)
        check(times.shape==(count+1,)
              and states.shape==(count+1,DOFS[mesh])
              and np.all(np.isfinite(states))
              and np.allclose(times,np.linspace(0,TF,count+1),
                              rtol=0,atol=5e-13),
              "Original native physical trajectory samples not on fixed time grid")
        traces={k:np.asarray(data[k],dtype=float) for k in required-{"times","states"}}
    for k,v in traces.items():
        check(v.shape==(count+1,) and np.all(np.isfinite(v)),
              f"Invalid original full physical trace {k}")
    # Redundancy: recompute directly from archived GetFEM energies and
    # original time-dependent h-weighted port/dissipation traces.
    h=traces["H"]
    supplied=np.r_[0.,np.cumsum(np.diff(times)*(
        traces["p_normal"][1:]+traces["p_viscous"][1:]))]
    dissipated=np.r_[0.,np.cumsum(np.diff(times)*(
        traces["d_grad"][1:]+traces["d_div"][1:]))]
    balance=h-h[0]+dissipated-supplied
    for name,real,reconstructed in (
        ("weighted wall work",traces["supplied"],supplied),
        ("viscous dissipation",traces["dissipated"],dissipated),
        ("physical Hamiltonian power",traces["balance"],balance),
    ):
        check(np.max(np.abs(real-reconstructed)) < 1e-12,
              f"Original native {name} differs from independent host integration")
    check(math.isclose(balance[-1],metrics["balance_defect"],
                       rel_tol=1e-7,abs_tol=2e-12),
          "Native archived h-weighted physical balance contradicts JSON")
    hfirst, hwidth = report["dynamic_variable_intervals"]["h"]
    pfirst, pwidth = report["dynamic_variable_intervals"]["p"]
    check(np.array_equal(
        states[0,hfirst:hfirst+hwidth],
        np.full(hwidth,5.,dtype=float)
    ) and np.array_equal(
        states[0,pfirst:pfirst+pwidth],
        np.zeros(pwidth,dtype=float)
    ),"Original prescribed native h0=5,p0=0 differs from archived physical state")
    ii=np.asarray([
        k for name in ("h","p")
        for k in range(
            report["dynamic_variable_intervals"][name][0],
            sum(report["dynamic_variable_intervals"][name])
        )
    ],dtype=np.int64)
    diff=norm(states[-1,ii]-states[0,ii])
    check(math.isclose(
        diff,float(report["final_dynamic_state_change_l2"]),
        rel_tol=1e-10,abs_tol=1e-12
    ),"Original dynamic physical state change differs from native FE trace")
    return {**metrics, "times":times, "states":states, "traces":traces,
            "final":states[-1], "initial":states[0],
            "final_dynamic":states[-1,ii].copy(), "dynamic_indices":ii}


def compare(directory):
    root=Path(directory)
    manifest=json.loads((root/"manifest.json").read_text(encoding="utf-8"))
    check(manifest.get("schema")=="v130-issue28-full-spinup-v1"
          and manifest.get("image_revision")==OCI
          and manifest.get("original_production_sources_unchanged") is True
          and manifest.get("fixed_physical_final_time")==TF
          and manifest.get("energy_convergence_claimed") is False,
          "Pinned native fixed-final-time refinement provenance mismatch")
    summary={}
    for mesh in MESHES:
        per_count={}
        for count in STEP_COUNTS:
            a=load_native_case(root,mesh,1,count)
            b=load_native_case(root,mesh,2,count)
            for key in ("Fa_max","TSCN_max","balance_defect",
                        "volume_defect","delta_H","W","D"):
                va,vb=a[key],b[key]
                tol=2e-12+5e-7*max(abs(va),abs(vb))
                check(abs(va-vb)<=tol,
                      f"Original native MPI physical diagnostic mismatch {key}")
            check(np.array_equal(a["initial"],b["initial"]),
                  "MPI1/2 published original analytic q0 differ bitwise")
            check(np.array_equal(
                a["dynamic_indices"],b["dynamic_indices"]
            ),"MPI1/2 original physical h/p intervals differ")
            dscale=norm(a["final"]-a["initial"])
            diff=norm(a["final"]-b["final"])
            check(diff<=1e-11+1e-8*dscale,
                  "Native MPI1/2 full final FE coefficient parity failed")
            per_count[count]=a
            print(
                f"mesh={mesh:g} n={count} dt={TF/count:.9g}: "
                f"TSCN_max={a['TSCN_max']:.9g} Fa_max={a['Fa_max']:.9g} "
                f"balance={a['balance_defect']:.9g} "
                f"balance/exchange={a['balance_over_exchange']:.9g} "
                f"volume={a['volume_defect']:.9g} MPI1/2 PASS",
                flush=True,
            )
        reference_indices=per_count[1]["dynamic_indices"]
        for count in STEP_COUNTS[1:]:
            check(np.array_equal(
                per_count[count]["dynamic_indices"],reference_indices
            ),"Native physical h/p FE intervals changed across temporal grids")
            check(np.array_equal(
                per_count[count]["initial"],per_count[1]["initial"]
            ),"Original prescribed initial physical state changed with dt")
        # Same spatial geometry and fixed final physical time.
        # Compare ONLY original dynamic h,p coefficients, not
        # algebraic multipliers that need not converge similarly.
        # These coefficient norms are not continuum L2 error norms.
        errors={}
        for low,high in zip(STEP_COUNTS[:-1],STEP_COUNTS[1:]):
            errors[f"{low}_vs_{high}"]=norm(
                per_count[low]["final_dynamic"]-per_count[high]["final_dynamic"]
            )
        orders={}
        vals=list(errors.values())
        for i in range(len(vals)-1):
            value=None
            if vals[i]>1e-13 and vals[i+1]>1e-13:
                value=math.log2(vals[i]/vals[i+1])
            orders[f"{STEP_COUNTS[i]}_{STEP_COUNTS[i+1]}_"
                   f"vs_{STEP_COUNTS[i+1]}_{STEP_COUNTS[i+2]}"]=value
        summary[f"h{mesh:g}"]={
            "final_dynamic_h_p_temporal_differences":errors,
            "observed_order_estimates_not_general_certificates":orders,
            "n10_balance_over_exchange":per_count[10]["balance_over_exchange"],
            "n40_balance_over_exchange":per_count[40]["balance_over_exchange"],
            "mpi_full_coefficient_parity":True,
            "endpoint_time_order_qualified":False,
        }
        print(f"mesh={mesh:g} temporal endpoint differences {errors}",flush=True)
        print(f"mesh={mesh:g} observed orders (diagnostic only) {orders}",flush=True)
    result={
        "schema":"scrimp-v130-issue28-full-spinup-compare-v1",
        "status":"NATIVE_TSCN_FULL_C1_RAMP_REFINEMENT_DIAGNOSTIC_ONLY",
        "cases":summary,
        "general_high_index_initializer_qualified":False,
        "nonlinear_mass_TSCN_general_accuracy_qualified":False,
        "multi_cycle_rotating_physical_SWE_qualified":False,
        "long_time_power_balance_qualified":False,
        "interpretation":(
            "Same physical final time and original nonlinear SWE on 2 meshes; "
            10,20,40 accepted PETSc TSCN steps independently verified "
            "against original per-step endpoint residuals and algebraic "
            "constraints. MPI1/2 full FE parity; true h-weighted work and "
            "viscous dissipation integrated independently. Endpoint "
            "refinement differences and order estimates are measured, "
            "not a blanket guarantee of nonlinear-mass TSCN accuracy."
        ),
    }
    (root/"full_spinup_comparison.json").write_text(
        json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print("COMPLETE: physical TSCN dt refinement diagnostic; NO general DAE/mass/power claim.")
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--evidence",type=Path,required=True)
    args=p.parse_args()
    compare(args.evidence)


if __name__=="__main__":
    main()
