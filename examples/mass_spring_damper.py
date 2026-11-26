# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Mass--spring--damper: a finite-dimensional port-Hamiltonian example.

State x=(q,p), with displacement q and momentum p=m*q_dot.
    H(q,p) = k*q**2/2 + p**2/(2*m)
    e_q = dH/dq = k*q,  e_p = dH/dp = p/m
    q_dot = e_p
    p_dot = -e_q - c*e_p + f(t)

Power balance: H_dot = -c*e_p**2 + f*e_p.

Each lumped port has a **single** DOF. The interval used here is only
a normalization domain for GetFEM's weak-form assembly; it is NOT a
spatial discretization of the oscillator.
"""

import numpy as np
import scrimp as S


def mass_spring_damper(
    final_time=8.0, time_step=0.01, plot=False, *,
    factor_solver=None, initialization="exact",
):
    """Simulate the forced, damped oscillator using lumped dynamic ports.

    Returns
    -------
    scrimp.DPHS
        Solved system with Hamiltonian terms and time histories.

    By default, impose the *exact* constitutive co-energies at t=0
    and preserve the prescribed differential state (q0,p0). The
    historical one-step pseudo-time relaxation is reproduced with
    initialization="legacy_pseudo"; it moves the physical initial
    momentum, and is therefore not the qualified default.
    """
    if initialization not in ("exact", "legacy_pseudo"):
        raise ValueError("initialization must be 'exact' or 'legacy_pseudo'")
    mass = 1.5
    stiffness = 5.0
    damping = 0.25
    interval_length = 1.0

    sys = S.DPHS("real")
    sys.set_domain(S.Domain("Interval", {"L": interval_length, "h": 1.0}))

    q = S.State("q", "Displacement", "lumped")
    p = S.State("p", "Linear momentum", "lumped")
    for state in (q, p):
        sys.add_state(state)
    sys.add_costate(S.CoState("e_q", "Spring force", q))
    sys.add_costate(S.CoState("e_p", "Velocity", p))

    # Unlike FEM("q", ...), this registers a single scalar DOF.
    sys.add_lumped_port("q")
    sys.add_lumped_port("p")

    for key, val in {
        "ell": interval_length,
        "osc_mass": mass,
        "spring_k": stiffness,
        "damp_c": damping,
    }.items():
        sys.gf_model.add_macro(key, str(float(val)))
    sys.gf_model.add_macro("force", "0.2*sin(2*pi*t)")

    bricks = [
        S.Brick("M_q", "q*Test_q/ell", [1], dt=True, position="flow"),
        S.Brick("M_p", "p*Test_p/ell", [1], dt=True, position="flow"),
        # Canonical Dirac interconnection J = [[0,1],[-1,0]]
        S.Brick("J_qp", "e_p*Test_q/ell", [1], position="effort"),
        S.Brick("J_pq", "-e_q*Test_p/ell", [1], position="effort"),
        S.Brick("Damping", "-damp_c*e_p*Test_p/ell",
                [1], position="effort"),
        # A time-only source is classified non-linear so that it is
        # assembled on every residual evaluation.
        S.Brick("External_force", "force*Test_p/ell", [1],
                linear=False, position="effort"),
        S.Brick("Spring_law", "(spring_k*q-e_q)*Test_e_q/ell", [1]),
        S.Brick("Momentum_law", "(p/osc_mass-e_p)*Test_e_p/ell", [1]),
    ]
    for brick in bricks:
        sys.add_brick(brick)

    q0, p0 = 0.4, 0.0
    sys.set_initial_value("q", q0)
    sys.set_initial_value("p", p0)
    if initialization == "exact":
        # These are true zero-dimensional GetFEM variables (one DOF
        # each), NOT FE fields to interpolate. The algebraic laws are
        # e_q = dH/dq = k*q and e_p = dH/dp = p/m. Initialize both
        # exactly rather than evolving (q,p) during a pseudo-TS.
        sys.gf_model.set_variable("e_q", np.array([stiffness*q0]))
        sys.gf_model.set_variable("e_p", np.array([p0/mass]))
    # Explicitly expose the mode for physical acceptance diagnostics.
    sys.mass_spring_initialization = initialization
    # An omitted factor_solver leaves the current SCRIMP/PETSc default
    # (MUMPS) untouched. Users can explicitly select superlu (serial),
    # superlu_dist (MPI), or another PETSc-supported LU backend.
    scheme = {}
    if factor_solver is not None:
        scheme["pc_factor_mat_solver_type"] = factor_solver
    # With fully consistent initial differential AND algebraic DOFs
    # no initial TS is necessary. Keep historical pseudo only by
    # explicit opt-in; the physical BDF2 time integrator is unchanged.
    init_options = (
        {"init_step": False}
        if initialization == "exact"
        else {"init_step_ts_type": "pseudo", "init_step": True}
    )
    sys.set_time_scheme(
        ts_type="bdf", ts_bdf_order=2,
        ts_exact_final_time="matchstep",
        ts_adapt_type="none",
        ksp_type="gmres", pc_type="lu",
        t_0=0.0, t_f=float(final_time),
        dt=float(time_step), dt_save=float(time_step),
        **init_options,
        **scheme,
    )
    sys.solve()

    sys.hamiltonian.set_name("Mass--spring--damper energy")
    sys.hamiltonian.add_term(
        S.Term("Spring potential", "spring_k*q*q/(2*ell)", [1])
    )
    sys.hamiltonian.add_term(
        S.Term("Kinetic", "p*p/(2*osc_mass*ell)", [1])
    )

    if plot:
        import matplotlib.pyplot as plt

        times = np.asarray(sys.solution["t"], dtype=float)
        displacement = np.asarray(sys.get_solution("q"), dtype=float).ravel()
        momentum = np.asarray(sys.get_solution("p"), dtype=float).ravel()
        fig, ax = plt.subplots()
        ax.plot(times, displacement, label="q (displacement)")
        ax.plot(times, momentum/mass, label="q_dot (velocity)")
        ax.set(xlabel="Time", title="Mass--spring--damper")
        ax.grid(True)
        ax.legend()
        fig.tight_layout()
        sys.plot_Hamiltonian(with_powers=False)
        plt.show()
    return sys


if __name__ == "__main__":
    mass_spring_damper(plot=True)
