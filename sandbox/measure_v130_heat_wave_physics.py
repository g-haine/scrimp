# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Non-invasive physical observables for the original 2D heat-wave example.

The saved accepted GetFEM states are evaluated with SCRIMP's own
get_quantity(), preserving the original mass/Dirac/gyrator weak forms.
Model verification, not a replacement coupling or a forced PASS.

For the corrected unit-parameter co-energy forms, the intended
semidiscrete closed-system identity is Hdot = -D_Q, where
H = 0.5*int_heat T^2 + 0.5*int_wave (p^2 + q.q),
D_Q = int_heat e_Q.e_Q >= 0,
and the two *physically signed* interface powers must cancel.
With the corrected temperature-input/traction-input weak forms:
P_interface_total=-orientation*U_T*Y_T+U_w*Y_w.
CN samples at each state mean are step powers; integrate dt*P[k]
rather than nodal trapezoids. This is NOT a BDF work identity.
"""

import numpy as np


REQUIRED = {
    "H_heat": (".5*T*T", 1),
    "H_wave_kinetic": (".5*p*p", 2),
    "H_wave_potential": (".5*q.q", 2),
    "D_heat_flux": ("e_Q.e_Q", 1),
    "P_interface_heat": ("U_T*Y_T", 10),
    "P_interface_wave": ("U_w*Y_w", 10),
    "P_outer": ("U_bnd*Y_bnd", 20),
}
# These surface trace diagnostics can be sensitive to FE traces in
# mixed/divergence-conforming spaces. Failure of a *diagnostic*
# should be recorded, not silently interpreted as zero.
OPTIONAL_TRACES = {
    "interface_U_T_minus_Y_w_L2_sq":
        ("(U_T-Y_w)*(U_T-Y_w)", 10),
    "interface_U_w_minus_Y_T_L2_sq":
        ("(U_w-Y_T)*(U_w-Y_T)", 10),
    "interface_U_w_plus_Y_T_L2_sq":
        ("(U_w+Y_T)*(U_w+Y_T)", 10),
    "interface_T_minus_p_L2_sq":
        ("(T-p)*(T-p)", 10),
    "interface_T_plus_p_L2_sq":
        ("(T+p)*(T+p)", 10),
    # With the common GetFEM interface Normal, the two outward
    # physical subdomain normals are opposite; flux continuity
    # corresponds to q.Normal + e_Q.Normal = 0 (not subtraction).
    "interface_normal_flux_sum_L2_sq":
        ("(q.Normal+e_Q.Normal)*(q.Normal+e_Q.Normal)",10),
}


def as_trace(values, times, name):
    arr = np.asarray(values)
    if np.iscomplexobj(arr):
        if np.max(abs(arr.imag)) > 1e-10:
            raise ValueError(f"{name}: nontrivial imaginary part")
        arr = arr.real
    arr = np.asarray(arr, dtype=float)
    if arr.shape != np.asarray(times).shape or not np.all(np.isfinite(arr)):
        raise ValueError(f"{name}: invalid original physical trace")
    return arr


def trapezoid(times, values):
    """Nodal work integral ONLY, never an asserted BDF2 stage identity."""
    t = np.asarray(times, dtype=float)
    p = as_trace(values, t, "quadrature")
    dt = np.diff(t)
    if t.ndim != 1 or t.size < 2 or not np.all(dt > 0.):
        raise ValueError("Original accepted-state time grid not increasing")
    return np.r_[0., np.cumsum(.5*dt*(p[1:]+p[:-1]))]


def collect_heat_wave_physics(
    system, times, H_total, *, heat_region=1, wave_region=2,
):
    """Actual original GetFEM/L2 powers and separate energy contributions.

    Never overwrites the original state, weak forms, or source examples.
    Optional interface *field* trace evaluations are recorded as errors
    if unsupported; algebraic port pairings remain required.
    """
    times = np.asarray(times, dtype=float)
    H_total = as_trace(H_total, times, "original total H")
    traces = {}
    if {heat_region,wave_region}!={1,2}:
        raise ValueError("Original Concentric model regions required")
    # Respect which of the two components occupies the inner disk.
    heat_orientation=1. if heat_region==1 else -1.
    dynamic_regions={
        "H_heat":heat_region, "H_wave_kinetic":wave_region,
        "H_wave_potential":wave_region, "D_heat_flux":heat_region,
    }
    for name, (expression, region) in REQUIRED.items():
        region=dynamic_regions.get(name,region)
        traces[name] = as_trace(
            system.get_quantity(expression, region=region),
            times, name
        )
    for name, (expression, region) in OPTIONAL_TRACES.items():
        try:
            traces[name] = as_trace(
                system.get_quantity(expression, region=region),
                times, name
            )
        except Exception as exc:
            traces[name] = None
            traces[name + "_error"] = (
                f"{type(exc).__name__}: {str(exc)[:220]}"
            )

    # These are the true GetFEM step powers of the production CN
    # scheme, evaluated at 0.5*(z[k-1]+z[k]) (and at t[k]).
    # The k=0 sample is instantaneous and is NOT a completed step.
    for name,(expression,region) in (
        ("CN_D_heat_flux",("e_Q.e_Q",heat_region)),
        ("CN_P_interface_heat",("U_T*Y_T",10)),
        ("CN_P_interface_wave",("U_w*Y_w",10)),
        ("CN_P_outer",("U_bnd*Y_bnd",20)),
    ):
        traces[name]=as_trace(
            system.get_quantity(expression,region=region,CN=True),
            times,name
        )
    H_components = (
        traces["H_heat"]+traces["H_wave_kinetic"]
        +traces["H_wave_potential"]
    )
    D=traces["D_heat_flux"]
    P1=traces["P_interface_heat"]
    P2=traces["P_interface_wave"]
    Pb=traces["P_outer"]
    # Explicit physical sign: the heat temperature input carries
    # work -U_T*J_Q.n_heat, whereas the wave traction input carries
    # work +U_w*p. GetFEM's interface Normal is common to the
    # two subdomains, so restore the actual heat outward normal.
    total_port=-heat_orientation*P1+P2+Pb
    work_heat=trapezoid(times,P1)
    work_wave=trapezoid(times,P2+Pb)
    dissipated=trapezoid(times,D)
    heat_defect=(traces["H_heat"]-traces["H_heat"][0]
                 +dissipated+heat_orientation*work_heat)
    H_wave=traces["H_wave_kinetic"]+traces["H_wave_potential"]
    wave_defect=H_wave-H_wave[0]-work_wave
    total_defect=H_total-H_total[0]+dissipated-trapezoid(times,total_port)

    dt=np.diff(times)
    step_D=traces["CN_D_heat_flux"][1:]
    step_P_heat=traces["CN_P_interface_heat"][1:]
    step_P_wave=traces["CN_P_interface_wave"][1:]
    step_P_outer=traces["CN_P_outer"][1:]
    stage_physical_exchange=(
        -heat_orientation*step_P_heat+step_P_wave+step_P_outer
    )
    # Actual partitioned CN identity (not a BDF endpoint approximation).
    stage_energy_defect=np.diff(H_total)+dt*(step_D-stage_physical_exchange)
    stage_heat_defect=np.diff(traces["H_heat"])+dt*(
        step_D+heat_orientation*step_P_heat
    )
    stage_wave_defect=np.diff(H_wave)-dt*(step_P_wave+step_P_outer)
    stage_work=np.r_[0.,np.cumsum(dt*stage_physical_exchange)]
    stage_dissipation=np.r_[0.,np.cumsum(dt*step_D)]
    stage_cumulative_residual=(
        H_total-H_total[0]+stage_dissipation-stage_work
    )

    traces.update({
        "H_original":H_total,
        "H_components":H_components,
        "H_wave":H_wave,
        "P_all_ports":total_port,
        "P_interface_sum":P2-heat_orientation*P1,
        "CN_P_interface_physical":
            -heat_orientation*traces["CN_P_interface_heat"]
            +traces["CN_P_interface_wave"],
        "CN_P_all_ports_physical":
            -heat_orientation*traces["CN_P_interface_heat"]
            +traces["CN_P_interface_wave"]+traces["CN_P_outer"],
        "CN_per_step_H_balance_error":np.r_[0.,stage_energy_defect],
        "CN_per_step_heat_balance_error":np.r_[0.,stage_heat_defect],
        "CN_per_step_wave_balance_error":np.r_[0.,stage_wave_defect],
        "CN_cumulative_total_balance_error":stage_cumulative_residual,
        "CN_cumulative_dissipation":stage_dissipation,
        "CN_cumulative_physical_boundary_work":stage_work,
        "cum_heat_dissipation_nodal_trapz":dissipated,
        "cum_heat_interface_work_nodal_trapz":work_heat,
        "cum_wave_interface_work_nodal_trapz":work_wave,
        "heat_balance_nodal_trapz_defect":heat_defect,
        "wave_balance_nodal_trapz_defect":wave_defect,
        "total_balance_nodal_trapz_defect":total_defect,
    })
    numerical={k:v for k,v in traces.items() if isinstance(v,np.ndarray)}
    extras={k:v for k,v in traces.items() if isinstance(v,str)}
    if np.min(D) < -1e-10:
        raise ValueError("Original GetFEM heat dissipation is negative")
    H_mismatch=float(np.max(abs(H_components-H_total)))
    if H_mismatch>1e-7*max(1.,float(np.max(abs(H_total)))):
        raise ValueError("Original Hamiltonian and GetFEM energy terms disagree")
    scale=max(1.,float(np.max(abs(P1))),float(np.max(abs(P2))))
    report={
        "schema":"scrimp-v130-original-heat-wave-physical-observables-v1",
        "scientific_status":"DIAGNOSTIC_ONLY_MODEL_NOT_YET_QUALIFIED",
        "source":"original_heat_wave_eq_unchanged",
        "Hamiltonian_terms_max_difference":H_mismatch,
        "dissipation_min":float(np.min(D)),
        "dissipation_max":float(np.max(D)),
        "max_interface_port_power_sum":
            float(np.max(abs(P2-heat_orientation*P1))),
        "relative_interface_port_power_sum":
            float(np.max(abs(P2-heat_orientation*P1))/scale),
        "max_CN_interface_power_defect":
            float(np.max(abs(
                traces["CN_P_interface_wave"]
                -heat_orientation*traces["CN_P_interface_heat"]
            ))),
        "max_CN_step_H_balance_defect":
            float(np.max(abs(stage_energy_defect))),
        "max_CN_step_heat_balance_defect":
            float(np.max(abs(stage_heat_defect))),
        "max_CN_step_wave_balance_defect":
            float(np.max(abs(stage_wave_defect))),
        "max_CN_cumulative_H_balance_defect":
            float(np.max(abs(stage_cumulative_residual))),
        "max_CN_positive_H_increment":
            float(np.max(np.diff(H_total))),
        "heat_normal_orientation":int(heat_orientation),
        "max_outer_port_power":float(np.max(abs(Pb))),
        "max_positive_total_H_increment":
            float(np.max(np.diff(H_total))),
        "max_abs_heat_balance_nodal_trapz_defect":
            float(np.max(abs(heat_defect))),
        "max_abs_wave_balance_nodal_trapz_defect":
            float(np.max(abs(wave_defect))),
        "max_abs_total_balance_nodal_trapz_defect":
            float(np.max(abs(total_defect))),
        "initial_H_total":float(H_total[0]),
        "final_H_total":float(H_total[-1]),
        "final_heat_dissipation_nodal_trapz":float(dissipated[-1]),
        "CN_step_powers_are_assembled_at_mean_FE_state":True,
        "nodal_trapezoid_is_diagnostic_only":True,
        "optional_surface_traces_unavailable":extras,
    }
    for key in OPTIONAL_TRACES:
        if isinstance(traces[key],np.ndarray):
            report[key+"_initial"]=float(traces[key][0])
            report[key+"_final"]=float(traces[key][-1])
    return numerical,report
