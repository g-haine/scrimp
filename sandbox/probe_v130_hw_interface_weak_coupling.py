# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native Heat-Wave interface port/amplitude/weak-coupling audit, MPI1.

Original published heat-wave, original GetFEM mixed weak forms,
strict algebraic DAE, corrected p exterior trace, real MUMPS CN and
original physical stage-mean energy observables are unchanged.

A deliberately narrow collector wrapper inspects the *same accepted*
solution for norms of port variables, RMS local power densities (to
distinguish near-zero field from cancellation of signed power), and
actual PETSc sparse stiffness blocks connecting physical states to
interface ports.

DO NOT infer strong trace continuity from raw DG/RT traces: the
physical constraints are weak and depend on chosen port causality.
This is exploratory native evidence, NOT a release model certificate.
"""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import csr_matrix

import probe_v130_main_2d_examples as base
import probe_v130_heat_wave_compatible_cn as physical_cn


N = 7689
REGION_INTERFACE = 10
BLOCKS = {
    "p_to_Y_w_trace": ("Y_w", "p"),
    "Y_w_to_U_T_feedback": ("U_T", "Y_w"),
    "e_Q_to_Y_T_flux": ("Y_T", "e_Q"),
    "Y_T_to_U_w_feedback": ("U_w", "Y_T"),
    "U_T_to_heat_flux_equation": ("e_Q", "U_T"),
    "U_w_to_wave_velocity_equation": ("p", "U_w"),
    "p_to_wave_stress_equation": ("q", "p"),
    "e_Q_to_heat_temperature_equation": ("T", "e_Q"),
}
FIELDS = (
    "U_T", "Y_w", "U_w", "Y_T", "T", "p",
)
FIELD_SQUARED = {
    name: f"{name}*{name}" for name in FIELDS
}
FIELD_SQUARED.update({
    "heat_normal_flux": "(e_Q.Normal)*(e_Q.Normal)",
    "wave_normal_stress": "(q.Normal)*(q.Normal)",
    "heat_interface_local_power_density_squared":
        "(U_T*Y_T)*(U_T*Y_T)",
    "wave_interface_local_power_density_squared":
        "(U_w*Y_w)*(U_w*Y_w)",
})


def sparse_coupling_audit(system):
    """Actual native assembled GetFEM/PETSc matrix coupling in FE order."""
    matrix = system.stiffness
    if matrix.getSize() != (N, N):
        raise RuntimeError("Original assembled Heat-Wave stiffness changed")
    indptr,indices,values=matrix.getValuesCSR()
    K=csr_matrix(
        (np.asarray(values).copy(),np.asarray(indices).copy(),
         np.asarray(indptr).copy()),
        shape=(N,N),
    )
    names=("T","p","q","e_Q","U_T","Y_T","U_w","Y_w",
           "U_bnd","Y_bnd")
    spans={}
    for name in names:
        lo,count=map(int,system.gf_model.interval_of_variable(name))
        if lo<0 or count<1 or lo+count>N:
            raise RuntimeError("Invalid GetFEM interval: "+name)
        spans[name]=[lo,count]
    full=[]
    for lo,n in spans.values():
        full.extend(range(lo,lo+n))
    if sorted(full)!=list(range(N)):
        raise RuntimeError("Port intervals overlap or omit original DOFs")
    stats={}
    for label,(row,col) in BLOCKS.items():
        rl,rn=spans[row]
        cl,cn=spans[col]
        block=K[rl:rl+rn,cl:cl+cn].tocsr()
        data=np.asarray(block.data)
        stats[label]={
            "row_variable":row,
            "column_variable":col,
            "shape":list(block.shape),
            "stored_entries":int(block.nnz),
            "numerically_nonzero_entries":int(np.count_nonzero(data)),
            "Frobenius_norm":float(np.linalg.norm(data)),
            "max_absolute_entry":float(np.max(np.abs(data))) if len(data) else 0.,
        }
    algebraic=("Y_w","Y_T","U_w","U_T","U_bnd","Y_bnd")
    raw=np.asarray(system.solution["z"])
    if raw.shape!=(5,N):
        raise RuntimeError("Expected original 4 CN steps, full GetFEM state")
    state_norms={}
    for name,(lo,count) in spans.items():
        coeff=raw[:,lo:lo+count]
        state_norms[name]={
            "coefficient_l2_by_saved_time":
                np.linalg.norm(coeff,axis=1).tolist(),
            "coefficient_max_absolute_by_saved_time":
                np.max(np.abs(coeff),axis=1).tolist(),
        }
    # All algebraic rows of the actual K z should vanish to
    # numerical accuracy for the homogeneous original model.
    # This is a diagnostic, not a replacement for actual IFunction.
    algebraic_weak_rows={}
    for name in algebraic:
        lo,count=spans[name]
        weak=np.asarray((K@raw.T)[lo:lo+count,:]).T
        algebraic_weak_rows[name]={
            "Kz_l2_by_saved_time":
                np.linalg.norm(weak,axis=1).tolist(),
            "Kz_max_absolute_by_saved_time":
                np.max(np.abs(weak),axis=1).tolist(),
        }
    return {
        "full_state_DOFS":N,
        "stiffness_stored_entries":int(K.nnz),
        "intervals":spans,
        "coupling_blocks":stats,
        "port_coefficient_norms":state_norms,
        "algebraic_Kz_row_residuals":algebraic_weak_rows,
        "not_a_strong_trace_constraint_proof":True,
    }


def surface_diagnostics(system,times):
    """True GetFEM region-10 quantities on original accepted CN states."""
    times=np.asarray(times)
    data={}
    failures={}
    for name,expr in FIELD_SQUARED.items():
        for midpoint in (False,True):
            label=name+("_CN" if midpoint else "_nodes")
            try:
                output=np.asarray(system.get_quantity(
                    expr,region=REGION_INTERFACE,CN=midpoint
                ),dtype=float).reshape(-1)
                if output.shape!=times.shape or np.any(~np.isfinite(output)):
                    raise RuntimeError("Nonfinite/incorrect GetFEM FE surface trace")
                # Quadratic functions must integrate nonnegatively
                # up to small absolute floating roundoff.
                if np.min(output)<-1e-9:
                    raise RuntimeError("Unexpected negative squared field integral")
                data[label]=output.tolist()
            except Exception as exc:
                failures[label]=f"{type(exc).__name__}: {str(exc)[:500]}"
    required=("U_T_CN","Y_w_CN","U_w_CN","Y_T_CN",
              "p_CN","T_CN")
    if any(name not in data for name in required):
        raise RuntimeError(
            "Required port/physical squared-norm GetFEM traces failed: "+
            json.dumps(failures,sort_keys=True)
        )
    # Power density squared is a separate, especially valuable
    # cancellation diagnostic, but may not be representable in
    # every GetFEM GWFL grammar. Do not report zero on failure.
    paired={}
    for pair in ("heat","wave"):
        key=f"{pair}_interface_local_power_density_squared_CN"
        squared=data.get(key)
        paired[pair]={
            "pointwise_power_squared_integral_by_CN_step":squared,
            "diagnostic_available":squared is not None,
        }
    return data, failures, paired


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--candidate",choices=("clamp","constrained_l2"),
                        default="clamp")
    args=parser.parse_args()
    if args.out.exists() and any(args.out.iterdir()):
        raise RuntimeError("New blank evidence directory required")
    original_collector=base.collect_heat_wave_physics
    captured={}
    invocation=0

    def observe_original_accepted_states(system,times,H,*a,**kw):
        nonlocal invocation
        invocation+=1
        if invocation!=1:
            raise RuntimeError("Expected exactly one original physical collector")
        # Preserve the original scientific observations and returned data.
        traces, report=original_collector(system,times,H,*a,**kw)
        matrix=sparse_coupling_audit(system)
        squared,failed,paired=surface_diagnostics(system,times)
        if not np.allclose(times,np.arange(5)*0.001,atol=1e-12,rtol=0):
            raise RuntimeError("Original CN time grid changed")
        # Recheck the weak power-pair equality in the archived
        # production observables before interpreting port amplitudes.
        signed_heat=np.asarray(traces["CN_P_interface_heat"])
        signed_wave=np.asarray(traces["CN_P_interface_wave"])
        captured.update({
            "schema":"scrimp-v130-original-heat-wave-interface-observability-v1",
            "status":"DIAGNOSTIC_ONLY_INTERFACE_COUPLING_PHYSICS_OPEN",
            "candidate":args.candidate,
            "source_example_unchanged":True,
            "strict_DAE_and_original_main_CN_executed":True,
            "original_physical_observables_unmodified":True,
            "full_4_step_CN_TIMES":np.asarray(times).tolist(),
            "original_max_signed_CN_interface_power":float(max(
                np.max(np.abs(signed_heat)),
                np.max(np.abs(signed_wave))
            )),
            "original_max_CN_power_exchange_defect":float(
                np.max(np.abs(signed_heat-signed_wave))
            ),
            "original_initial_T_minus_p_trace_L2_sq":
                report["interface_T_minus_p_L2_sq_initial"],
            "original_initial_flux_sum_trace_L2_sq":
                report["interface_normal_flux_sum_L2_sq_initial"],
            "matrix":matrix,
            "true_region10_squared_integrals":squared,
            "optional_GWFL_trace_errors":failed,
            "power_density_squared":paired,
            "raw_trace_mismatch_is_not_proof_of_weak_FEM_constraint_failure":True,
            "no_published_example_or_production_code_change":True,
            "full_multiphysics_coupling_qualified":False,
            "MPI2_qualified":False,
        })
        file=args.out/"interface_weak_coupling_diagnostic.json"
        if file.exists():
            raise RuntimeError("Original interface audit must not overwrite evidence")
        file.write_text(json.dumps(captured,indent=2)+"\n",encoding="utf-8")
        print("INTERFACE AUDIT saved: "+str(file),flush=True)
        for label,info in matrix["coupling_blocks"].items():
            print(
                "  MATRIX "+label+
                f" nonzero={info['numerically_nonzero_entries']} "+
                f"frob={info['Frobenius_norm']:.8g}",
                flush=True
            )
        for name in ("U_T_CN","Y_w_CN","U_w_CN","Y_T_CN","p_CN","T_CN",
                     "heat_interface_local_power_density_squared_CN",
                     "wave_interface_local_power_density_squared_CN"):
            if name in squared:
                print(
                    f"  SURFACE {name}={squared[name]}",
                    flush=True,
                )
            else:
                print(f"  SURFACE {name}: UNAVAILABLE "+
                      str(failed.get(name)),flush=True)
        return traces,report

    old_argv=sys.argv
    base.collect_heat_wave_physics=observe_original_accepted_states
    try:
        sys.argv=[
            str(Path(__file__).name),
            "--candidate",args.candidate,
            "--backend","mumps",
            "--out",str(args.out),
        ]
        physical_cn.main()
    finally:
        sys.argv=old_argv
        base.collect_heat_wave_physics=original_collector
    if invocation!=1 or not captured:
        raise RuntimeError("Actual original physical trace collector not executed")
    print(
        "PASS_INTERFACE_DIAGNOSTIC_COLLECTED_COUPLING_PHYSICS_NOT_QUALIFIED",
        flush=True,
    )


if __name__=="__main__":
    main()
