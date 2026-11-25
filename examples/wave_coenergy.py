# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             examples/wave_coenergy.py
- authors:          Giuseppe Ferraro, Ghislain Haine
- date:             22 nov. 2022
- brief:            wave equations in co-energy formulation, two sub-domains
"""

# Import scrimp
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import scrimp as S

# TUTORIAL 2C-00: same T and rho as 2A/2B; true inverse compliance
# GetFEM nested tensors encode columns. T is symmetric, so its inverse is
# [[T22/det, -T12/det], [-T12/det, T11/det]].
DET_T = "((5+x)*(2+y)-(x*y)*(x*y))"
TINV = (
    "[[(2+y)/" + DET_T + ",-(x*y)/" + DET_T + "],"
    "[-(x*y)/" + DET_T + ",(5+x)/" + DET_T + "]]"
)
FEEDBACK_GAIN = 0.5


def wave_coenergy_eq(
    *, mesh_size=0.2, final_time=0.08, time_step=0.008,
    factor_solver="mumps",
):
    """Solve the anisotropic 2D wave in genuine co-energy variables.

    q=e_q (stress), p=e_p (velocity); Tinv=T**(-1), rho as in 2A/2B.
    Damping acts only on the inner disk, feedback only on the outer rim.
    The short fixed CN scenario is for native acceptance, NOT convergence.
    """
    mesh_size = float(mesh_size)
    final_time = float(final_time)
    time_step = float(time_step)
    if not (np.isfinite(mesh_size) and mesh_size>0
            and np.isfinite(final_time) and final_time>0
            and np.isfinite(time_step) and time_step>0):
        raise ValueError("Mesh and time controls must be positive and finite")
    ratio = final_time/time_step
    if abs(ratio-round(ratio))>1e-10:
        raise ValueError("Use an integral number of fixed CN steps")
    if factor_solver not in (None, "mumps", "superlu", "superlu_dist"):
        raise ValueError("Unknown PETSc LU solver")


    # Init the distributed port-Hamiltonian system
    wave = S.DPHS("real")

    # TUTORIAL 2C-01: concentric disk/annulus and physical labels
    # Set the domain (using the built-in geometry `Concentric`)
    # Labels: Disk = 1, Annulus = 2, Interface = 10, Boundary = 20
    omega = S.Domain("Concentric", {"R": 1.0, "r": 0.6, "h": mesh_size})
    
    # And add it to the dphs
    wave.set_domain(omega)

    # TUTORIAL 2C-02: substituted coenergy states
    ## Define the variables
    states = [
        S.State("q", "Stress", "vector-field"),
        S.State("p", "Velocity", "scalar-field"),
    ]
    # Use of the `substituted=True` keyword to get the co-energy formulation
    costates = [
        S.CoState("e_q", "Stress", states[0], substituted=True),
        S.CoState("e_p", "Velocity", states[1], substituted=True),
    ]
    
    # Add them to the dphs
    for state in states:
        wave.add_state(state)
    for costate in costates:
        wave.add_costate(costate)
    
    # TUTORIAL 2C-03: inner resistive port and outer boundary port
    # Define the dissipative port, only on the subdomain labelled 1 = the internal disk
    ports = [
        S.Port("Damping", "e_r", "e_r", "scalar-field", substituted=True, region=1),
    ]
    
    # Add it to the dphs
    for port in ports:
        wave.add_port(port)
    
    # Define the control port
    control_ports = [
        S.Control_Port(
            "Boundary control",
            "U",
            "Normal force",
            "Y",
            "Velocity trace",
            "scalar-field",
            region=20,
        ),
    ]

    # Add it to the dphs
    for ctrl_port in control_ports:
        wave.add_control_port(ctrl_port)
    
    # Define the Finite Elements Method of each port
    FEMs = [
        S.FEM(states[0].get_name(), 1, "DG"),
        S.FEM(states[1].get_name(), 2, "CG"),
        S.FEM(ports[0].get_name(), 1, "DG"),
        S.FEM(control_ports[0].get_name(), 1, "DG"),
    ]
    
    # Add them to the dphs
    for FEM in FEMs:
        wave.add_FEM(FEM)
    
    # TUTORIAL 2C-04: actual T inverse, density and bounded 1/nu
    # Define physical parameters: care must be taken, 
    # in the co-energy formulation, some parameters are
    # inverted in comparison to the classical formulation
    parameters = [
        S.Parameter(
            "Tinv",
            "Anisotropic compliance T inverse",
            "tensor-field",
            TINV,
            "q",
        ),
        S.Parameter("rho", "Mass density", "scalar-field", "3-x", "p"),
        S.Parameter(
            "nu",
            "Viscosity",
            "scalar-field",
            # Strictly positive on region 1 for the 1/nu mass brick.
            # The damping port does not exist on region 2.
            "0.2+10*(0.36-(x*x+y*y))",
            ports[0].get_name(),
        ),
    ]

    # Add them to the dphs
    for parameter in parameters:
        wave.add_parameter(parameter)

    # TUTORIAL 2C-05: mass tensors, Dirac interconnection, resistive port
    # Define the pHs via `Brick` == non-zero block matrices == variational terms
    # Since we use co-energy formulation, constitutive relations are already taken into
    # account in the mass matrices M_q and M_p
    bricks = [
        ## Define the Dirac structure
        # Define the mass matrices from the left-hand side: the `flow` part of the Dirac structure
        S.Brick("M_q", "q.Tinv.Test_q", [1, 2], dt=True, position="flow"),
        S.Brick("M_p", "p*rho*Test_p", [1, 2], dt=True, position="flow"),
        S.Brick("M_r", "e_r/nu*Test_e_r", [1], position="flow"),
        S.Brick("M_Y", "Y*Test_Y", [20], position="flow"),
        # Define the matrices from the right-hand side: the `effort` part of the Dirac structure
        S.Brick("D", "Grad(p).Test_q", [1, 2], position="effort"),
        S.Brick("-D^T", "-q.Grad(Test_p)", [1, 2], position="effort"),
        S.Brick("I_r", "e_r*Test_p", [1], position="effort"),
        S.Brick("B", "U*Test_p", [20], position="effort"),
        S.Brick("-I_r^T", "-p*Test_e_r", [1], position="effort"),
        S.Brick("-B^T", "-p*Test_Y", [20], position="effort"),
        ## Define the constitutive relations
        # Already taken into account in the Dirac Structure!
    ]
    
    # Add all these `Bricks` to the dphs
    for brick in bricks:
        wave.add_brick(brick)

    # TUTORIAL 2C-06: passive outer feedback and two true initial pulses
    ## Initialize the problem
    # The controls expression
    # U=-kY, Y=p: negative feedback power U*Y=-k*Y**2 <= 0.
    expressions = [f"-{FEEDBACK_GAIN}*Y"]
    
    # Add each expression to its control_port
    for control_port, expression in zip(control_ports, expressions):
        # Set the control functions (automatic construction of bricks such that -M_u u + f(t) = 0)
        wave.set_control(control_port.get_name(), expression)

    # Set the initial data
    wave.set_initial_value("q", "[0., 0.]")
    # Smooth two-part velocity pulse excites both internal damping and the
    # exterior feedback already at t=0, without prescribing a boundary trace.
    wave.set_initial_value(
        "p",
        "exp(-10*(x*x+y*y))+0.4*exp(-15*((x-0.8)*(x-0.8)+y*y))",
    )

    # TUTORIAL 2C-07: strict DAE initialization, fixed CN
    ## Solve in time
    # Define the time scheme ("cn" is Crank-Nicolson)
    options = {}
    if factor_solver is not None:
        options["pc_factor_mat_solver_type"] = factor_solver
    wave.set_time_scheme(
        init_step_ts_type="algebraic", init_step=True,
        ts_type="cn", ts_adapt_type="none",
        ts_exact_final_time="matchstep",
        ksp_type="gmres", pc_type="lu",
        t_0=0., t_f=final_time,
        dt=time_step, dt_save=time_step,
        **options,
    )

    # Solve
    wave.solve()

    # TUTORIAL 2C-08: genuine GetFEM coenergy Hamiltonian
    ## Post-processing
    ## Set Hamiltonian's name
    wave.hamiltonian.set_name("Mechanical energy")
    # Define each Hamiltonian Term
    terms = [
        S.Term("Potential energy", "0.5*q.Tinv.q", [1, 2]),
        S.Term("Kinetic energy", "0.5*p*p*rho", [1, 2]),
    ]
    # Add them to the Hamiltonian
    for term in terms:
        wave.hamiltonian.add_term(term)

    return wave


# TUTORIAL 2C-09: actual CN-midpoint FE work and dissipation audit
def measure_coenergy(system):
    """Audit actual CN midpoint FE work, damping and boundary feedback.

    The substituted resistive port has *no* independent flow unknown and
    Port.get_power() correctly cannot report its dissipative power.
    We assemble its true mixed FE dissipation in two mathematically
    equivalent ways: D_r=integral e_r**2/nu = -integral p*e_r.
    The pointwise nu*p**2 is ONLY a diagnostic because e_r is DG1 and
    nu*p is not generally in the same FE space.
    """
    import getfem as gf

    t = np.asarray(system.solution["t"], dtype=float)
    z = np.asarray(system.solution["z"], dtype=float)
    H = np.asarray(system.get_Hamiltonian(), dtype=float)
    if not system.powers_computed:
        system.compute_powers()
    boundary = system.ports["Boundary control"]
    P = np.asarray(boundary.get_power(), dtype=float)
    if boundary.get_power_sampling() != "cn_step":
        raise RuntimeError("True CN-step boundary power required")
    if (t.ndim!=1 or len(t)<2 or z.ndim!=2 or len(z)!=len(t)
            or H.shape!=t.shape or P.shape!=t.shape
            or not np.all(np.diff(t)>0)
            or any(not np.all(np.isfinite(a)) for a in (t,z,H,P))):
        raise RuntimeError("Invalid native coenergy CN FE trajectory")

    resistive = system.ports["Damping"]
    dr = np.empty(len(t))
    dr_cross = np.empty(len(t))
    dr_pointwise = np.empty(len(t))
    feedback = np.empty(len(t))
    feedback_constraint = np.empty(len(t))
    resistive_weak_residual_L2 = np.empty(len(t))

    def assemble(form, region, mesh_id):
        return float(gf.asm(
            "generic", system.domain._int_method[mesh_id],
            0, form, region, system.gf_model,
        ))

    for k, tk in enumerate(t):
        system.gf_model.set_time(float(tk))
        system.gf_model.to_variables(0.5*(z[k-1]+z[k]) if k else z[k])
        dr[k] = assemble("e_r*e_r/nu", 1, resistive.get_mesh_id())
        dr_cross[k] = assemble("-p*e_r", 1, resistive.get_mesh_id())
        dr_pointwise[k] = assemble("nu*p*p", 1, resistive.get_mesh_id())
        feedback[k] = assemble(
            f"{FEEDBACK_GAIN}*Y*Y", 20, boundary.get_mesh_id(),
        )
        feedback_constraint[k] = assemble(
            f"(U+{FEEDBACK_GAIN}*Y)*(U+{FEEDBACK_GAIN}*Y)",
            20, boundary.get_mesh_id(),
        )
        resistive_weak_residual_L2[k] = assemble(
            "(e_r/nu+p)*(e_r/nu+p)", 1, resistive.get_mesh_id(),
        )
    if any(not np.all(np.isfinite(a)) for a in
           (dr, dr_cross, dr_pointwise, feedback, feedback_constraint,
            resistive_weak_residual_L2)):
        raise RuntimeError("Nonfinite independent coenergy GetFEM quadrature")

    dt = np.diff(t)
    step_internal = dt*dr[1:]
    step_boundary = -dt*P[1:]
    Li = np.r_[0., np.cumsum(step_internal)]
    Lb = np.r_[0., np.cumsum(step_boundary)]
    balance = H-H[0]+Li+Lb
    defects = np.diff(H)+step_internal+step_boundary
    magnitude = max(1., float(np.max(np.abs(H))))
    init = system.initialization_report
    cv = (system.last_ts_converged_reason>0
          and system.last_ts_steps==len(t)-1
          and system.last_ts_rejections==0
          and system.last_ts_snes_failures==0)
    gates = {
        "native_fixed_CN_all_steps": bool(cv),
        "algebraic_DAE_initialization": bool(
            isinstance(init,dict) and init.get("status")=="PASS"
            and init.get("differential_initial_exact")),
        "true_compliance_positive_energy": bool(
            np.min(H)>-1e-9 and H[0]>1e-6),
        "strictly_positive_internal_dissipation": bool(
            np.min(dr[1:])>=-1e-10 and Li[-1]>1e-8),
        "negative_nontrivial_feedback_power": bool(
            np.max(P[1:])<=1e-9 and Lb[-1]>1e-8),
        "independent_weak_internal_loss_matches": bool(
            np.max(np.abs(dr-dr_cross))<=1e-8*magnitude),
        "independent_feedback_quadrature_matches": bool(
            np.max(np.abs(feedback+P))<=1e-8*magnitude),
        "feedback_algebraic_constraint": bool(
            np.max(np.abs(feedback_constraint))<=1e-8),
        "each_CN_step_dissipative_energy_balance": bool(
            np.max(np.abs(defects))<=1e-8*magnitude),
        "cumulative_CN_dissipative_energy_balance": bool(
            np.max(np.abs(balance))<=1e-8*magnitude),
    }
    return dict(
        t=t,z=z,H=H,power_boundary=P,
        dissipation_internal=dr,
        dissipation_internal_weak_check=dr_cross,
        dissipation_internal_pointwise_diagnostic=dr_pointwise,
        dissipation_feedback=feedback,
        feedback_constraint=feedback_constraint,
        resistive_weak_residual_L2=resistive_weak_residual_L2,
        cumulative_internal_loss=Li,cumulative_feedback_loss=Lb,
        energy_step_defects=defects,energy_cumulative_defects=balance,
        gates=gates,
    )


# TUTORIAL 2C-10: original FE trajectories, VTK and metadata
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver",default="mumps",
                        choices=("mumps","superlu","superlu_dist"))
    parser.add_argument("--mesh-size",type=float,default=.2)
    parser.add_argument("--final-time",type=float,default=.08)
    parser.add_argument("--time-step",type=float,default=.008)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--export",action="store_true")
    args = parser.parse_args()
    if (args.check or args.export) and args.output is None:
        parser.error("--check and --export require --output")
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    if comm.size>1 and args.solver=="superlu":
        parser.error("Sequential SuperLU cannot run in MPI2")
    system = wave_coenergy_eq(
        mesh_size=args.mesh_size,final_time=args.final_time,
        time_step=args.time_step,factor_solver=args.solver,
    )
    if args.output is None:
        return
    args.output.mkdir(parents=True,exist_ok=True)
    try:
        obs = measure_coenergy(system)
    except Exception as exc:
        np.savez_compressed(
            args.output/f"rank{comm.rank}.failed.npz",
            t=np.asarray(system.solution["t"]),
            z=np.asarray(system.solution["z"]),
        )
        (args.output/f"rank{comm.rank}.failed.json").write_text(
            json.dumps({"status":"UNQUALIFIED_NATIVE_COENERGY",
                        "error":f"{type(exc).__name__}: {exc}"},indent=2)+"\n"
        )
        raise
    np.savez_compressed(
        args.output/f"rank{comm.rank}.npz",
        **{k:v for k,v in obs.items() if isinstance(v,np.ndarray)},
    )
    meta = {
        "status":("PASS_NATIVE_WAVE2D_COENERGY_PHYSICS"
                  if all(obs["gates"].values())
                  else "FAIL_NATIVE_WAVE2D_COENERGY_PHYSICS"),
        "rank":comm.rank,"MPI":comm.size,
        "requested_LU":args.solver,
        "selected_LU":system.time_scheme.getString(
            "pc_factor_mat_solver_type","").lower(),
        "DOFs":int(obs["z"].shape[1]),
        "steps":int(system.last_ts_steps),
        "TS_reason":int(system.last_ts_converged_reason),
        "DAE_initialization":system.initialization_report,
        "physical_gates":obs["gates"],
        "max_abs_step_defect":float(np.max(np.abs(obs["energy_step_defects"]))),
        "max_abs_cumulative_defect":float(
            np.max(np.abs(obs["energy_cumulative_defects"]))),
        "total_internal_loss":float(obs["cumulative_internal_loss"][-1]),
        "total_feedback_loss":float(obs["cumulative_feedback_loss"][-1]),
        "state_SHA256":hashlib.sha256(
            np.ascontiguousarray(obs["z"]).tobytes()).hexdigest(),
    }
    (args.output/f"rank{comm.rank}.json").write_text(
        json.dumps(meta,indent=2,default=str)+"\n",encoding="utf-8",
    )
    if args.check and not all(obs["gates"].values()):
        raise RuntimeError("Coenergy physical acceptance FAILED: "+json.dumps(meta))
    if args.export:
        # Collective GetFEM export on *all* MPI ranks.
        for name in ("q","p","e_r"):
            system.export_to_pv(name,path=str(args.output),t="All")
    if comm.rank==0:
        print(json.dumps(meta,indent=2,default=str),flush=True)


if __name__ == "__main__":
    main()
