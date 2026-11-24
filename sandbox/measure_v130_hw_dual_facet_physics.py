# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Actual one-sided Heat-Wave CN midpoint port-power/energy diagnostic.

The corrected interface uses 2D oriented DISK faces 30 for heat
boundary effort/flow, and ANNULUS inner faces 31 for wave velocity
and stress control. The old doubled region 10 is deliberately NEVER
used for computing interface power.

This is NOT a physical certification: energy identities and interface
excitation must be independently validated on the corrected native
CN trajectories, and can fail. We preserve all FE state and backend
evidence without replacing the original weak forms.
"""

import numpy as np


HEAT_FACES = 30
WAVE_FACES = 31


def numeric(system, times, expr, region, *, cn=False):
    arr=np.asarray(
        system.get_quantity(expr,region=region,CN=cn),
        dtype=float,
    ).reshape(-1)
    if arr.shape!=times.shape or not np.all(np.isfinite(arr)):
        raise RuntimeError(
            "Invalid native GetFEM physical trace "
            f"{expr} on region {region}, CN={cn}"
        )
    return arr


def collect_heat_wave_dual_facet_physics(
    system, times, H_total, *, heat_faces=HEAT_FACES, wave_faces=WAVE_FACES,
):
    """Real GetFEM nodal energies and CN-midpoint separately sided powers.

    Optional explicit region IDs allow the experimental generic compiler
    to supply actual oriented faces. Historical callers retain the 30/31
    default unchanged; no mutation of module globals is necessary.
    """
    heat_faces = int(heat_faces)
    wave_faces = int(wave_faces)
    if heat_faces == wave_faces or 10 in (heat_faces, wave_faces):
        raise RuntimeError("Physical power requires two distinct one-sided regions")
    times=np.asarray(times,dtype=float).reshape(-1)
    H=np.asarray(H_total,dtype=float).reshape(-1)
    if (times.shape!=(5,) or H.shape!=times.shape
            or not np.allclose(times,np.arange(5)*.001,atol=1e-12,rtol=0)
            or not np.all(np.isfinite(H))):
        raise RuntimeError("Wrong original 4-step Heat-Wave physical history")
    available = {int(r) for r in np.asarray(
        system.domain.get_mesh()[0].regions(), dtype=np.int64
    ).reshape(-1)}
    if heat_faces not in available or wave_faces not in available:
        raise RuntimeError(
            "Oriented physical interface facet regions missing: "
            f"heat={heat_faces}, wave={wave_faces}, available={sorted(available)}"
        )
    definitions={
        "H_heat":(".5*T*T",1),
        "H_wave_kinetic":(".5*p*p",2),
        "H_wave_potential":(".5*q.q",2),
        "D_heat_flux":("e_Q.e_Q",1),
        "P_interface_heat":("U_T*Y_T",heat_faces),
        "P_interface_wave":("U_w*Y_w",wave_faces),
        "P_outer":("U_bnd*Y_bnd",20),
    }
    traces={
        name:numeric(system,times,expr,reg)
        for name,(expr,reg) in definitions.items()
    }
    for name,expr,reg in (
        ("CN_D_heat_flux","e_Q.e_Q",1),
        ("CN_P_interface_heat","U_T*Y_T",heat_faces),
        ("CN_P_interface_wave","U_w*Y_w",wave_faces),
        ("CN_P_outer","U_bnd*Y_bnd",20),
        ("CN_interface_heat_effort_squared","U_T*U_T",heat_faces),
        ("CN_interface_heat_output_squared","Y_T*Y_T",heat_faces),
        ("CN_interface_wave_effort_squared","U_w*U_w",wave_faces),
        ("CN_interface_wave_output_squared","Y_w*Y_w",wave_faces),
    ):
        traces[name]=numeric(system,times,expr,reg,cn=True)

    H_wave=traces["H_wave_kinetic"]+traces["H_wave_potential"]
    H_parts=traces["H_heat"]+H_wave
    dt=np.diff(times)
    # Original heat-region=1 disk orientation is its outward normal on 30.
    # Original model control interconnection has U_T=Y_w, U_w=Y_T.
    # The port-Hamiltonian power orientation follows the existing
    # model's weak-form convention: -U_T Y_T + U_w Y_w + P_outer.
    signed_heat=traces["CN_P_interface_heat"][1:]
    signed_wave=traces["CN_P_interface_wave"][1:]
    outer=traces["CN_P_outer"][1:]
    dissipation=traces["CN_D_heat_flux"][1:]
    port_exchange=-signed_heat+signed_wave+outer
    defect_heat=np.diff(traces["H_heat"])+dt*(dissipation+signed_heat)
    defect_wave=np.diff(H_wave)-dt*(signed_wave+outer)
    defect_total=np.diff(H)+dt*(dissipation-port_exchange)
    cumulative_dissipation=np.r_[0.,np.cumsum(dt*dissipation)]
    cumulative_ports=np.r_[0.,np.cumsum(dt*port_exchange)]
    trace_error=float(np.max(np.abs(H_parts-H)))
    pheat=np.abs(signed_heat)
    pwave=np.abs(signed_wave)
    max_power=max(float(np.max(pheat)),float(np.max(pwave)))
    max_diss=float(np.max(np.abs(dissipation)))
    ratio=max_power/max(max_diss,1e-30)
    # Detect port constraint differences in the two distinct oriented
    # FE surface regions. These expressions are diagnostic; their
    # optional failure never becomes a fabricated zero.
    optional={}
    unavailable={}
    for name,expr,region in (
        ("CN_U_T_minus_Y_w_L2_sq",
         "(U_T-Y_w)*(U_T-Y_w)",heat_faces),
        ("CN_U_w_minus_Y_T_L2_sq",
         "(U_w-Y_T)*(U_w-Y_T)",wave_faces),
        ("CN_heat_surface_power_density_squared",
         "(U_T*Y_T)*(U_T*Y_T)",heat_faces),
        ("CN_wave_surface_power_density_squared",
         "(U_w*Y_w)*(U_w*Y_w)",wave_faces),
    ):
        try:
            optional[name]=numeric(system,times,expr,region,cn=True)
        except Exception as exc:
            unavailable[name]=f"{type(exc).__name__}: {str(exc)[:250]}"
    traces.update(optional)
    traces.update({
        "H_original":H,
        "H_components":H_parts,
        "H_wave":H_wave,
        "CN_P_interface_physical":traces["CN_P_interface_wave"]-
                                   traces["CN_P_interface_heat"],
        "CN_P_all_ports_physical":traces["CN_P_interface_wave"]-
                                   traces["CN_P_interface_heat"]+
                                   traces["CN_P_outer"],
        "CN_per_step_H_balance_error":np.r_[0.,defect_total],
        "CN_per_step_heat_balance_error":np.r_[0.,defect_heat],
        "CN_per_step_wave_balance_error":np.r_[0.,defect_wave],
        "CN_cumulative_total_balance_error":
            H-H[0]+cumulative_dissipation-cumulative_ports,
        "CN_cumulative_dissipation":cumulative_dissipation,
        "CN_cumulative_physical_boundary_work":cumulative_ports,
    })
    if float(np.min(traces["D_heat_flux"])) < -1e-10:
        raise RuntimeError("Negative GetFEM heat dissipation")
    if trace_error>1e-7*max(1.,float(np.max(np.abs(H)))):
        raise RuntimeError("Heat + wave GetFEM H differs from Hamiltonian")
    result={
        "schema":"scrimp-v130-HW-dual-facet-CN-physical-audit-v1",
        "scientific_status":"DIAGNOSTIC_ONLY_CORRECTED_INTERFACE_NOT_YET_QUALIFIED",
        "interface_heat_region":heat_faces,
        "interface_wave_region":wave_faces,
        "old_doubled_interface_region_excluded_from_physical_power":True,
        "source":"unmodified original heat_wave_eq, runtime dual-face remapping",
        "Hamiltonian_terms_max_difference":trace_error,
        "dissipation_min":float(np.min(traces["D_heat_flux"])),
        "dissipation_max":float(np.max(traces["D_heat_flux"])),
        "max_interface_port_power_sum":
            float(np.max(np.abs(
                traces["P_interface_wave"]-traces["P_interface_heat"]
            ))),
        "max_CN_interface_power_defect":
            float(np.max(np.abs(signed_wave-signed_heat))),
        "max_CN_step_H_balance_defect":
            float(np.max(np.abs(defect_total))),
        "max_CN_step_heat_balance_defect":
            float(np.max(np.abs(defect_heat))),
        "max_CN_step_wave_balance_defect":
            float(np.max(np.abs(defect_wave))),
        "max_CN_cumulative_H_balance_defect":
            float(np.max(np.abs(traces[
                "CN_cumulative_total_balance_error"
            ]))),
        "max_CN_positive_H_increment":float(np.max(np.diff(H))),
        "max_outer_port_power":float(np.max(np.abs(traces["P_outer"]))),
        "max_abs_total_balance_nodal_trapz_defect":None,
        "CN_step_powers_are_assembled_at_mean_FE_state":True,
        "nodal_trapezoid_is_diagnostic_only":True,
        "max_individual_CN_interface_power":max_power,
        "max_CN_heat_dissipation":max_diss,
        "interface_exchange_ratio_to_heat_dissipation":ratio,
        "interface_exchange_nontrivial":bool(ratio>1e-8),
        "CN_interface_power_balance_weak_form_expected_sign":
            f"-heat_region{heat_faces} + wave_region{wave_faces}",
        "individual_CN_interface_power_heat":signed_heat.tolist(),
        "individual_CN_interface_power_wave":signed_wave.tolist(),
        "max_interface_heat_power_density_integral_l2_sq":(
            None if "CN_heat_surface_power_density_squared" not in traces
            else float(np.max(traces[
                "CN_heat_surface_power_density_squared"
            ]))
        ),
        "max_interface_wave_power_density_integral_l2_sq":(
            None if "CN_wave_surface_power_density_squared" not in traces
            else float(np.max(traces[
                "CN_wave_surface_power_density_squared"
            ]))
        ),
        "optional_surface_traces_unavailable":unavailable,
        "energy_physics_qualified":False,
        "MPI2_qualified":False,
    }
    print(
        "CORRECTED ONE-SIDED PORTS: "
        f"max_heat={float(np.max(pheat)):.7g} "
        f"max_wave={float(np.max(pwave)):.7g} "
        f"ratio_to_dissipation={ratio:.7g} "
        f"max_signed_exchange_defect={result['max_CN_interface_power_defect']:.7g}",
        flush=True,
    )
    print(
        "CORRECTED CN: "
        f"heat_defect={result['max_CN_step_heat_balance_defect']:.7g} "
        f"wave_defect={result['max_CN_step_wave_balance_defect']:.7g} "
        f"total_defect={result['max_CN_step_H_balance_defect']:.7g}",
        flush=True,
    )
    return traces,result
