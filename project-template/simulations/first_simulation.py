"""A short, genuine SCRIMP port-Hamiltonian simulation.

Mass–spring–damper with a single spatially lumped port for each state.
Based on the maintained SCRIMP example `examples/mass_spring_damper.py`.
Run from the workspace root:
    python simulations/first_simulation.py

Outputs: `outputs/oscillator.csv` on YOUR computer (not in the image).
This introductory model is serial. Use a distributed PDE example to
validate real MPI numerical scaling; `mpi_check.py` checks MPI startup.
"""

from pathlib import Path

import numpy as np
import scrimp as S
from petsc4py import PETSc


def main() -> None:
    if PETSc.COMM_WORLD.getSize() != 1:
        raise RuntimeError("This first oscillator is serial; use mpi_check.py for MPI startup.")

    system = S.DPHS("real")
    system.set_domain(S.Domain("Interval", {"L": 1.0, "h": 1.0}))

    q = S.State("q", "Displacement", "lumped")
    p = S.State("p", "Momentum", "lumped")
    for state in (q, p):
        system.add_state(state)
    system.add_costate(S.CoState("e_q", "Spring force", q))
    system.add_costate(S.CoState("e_p", "Velocity", p))
    system.add_lumped_port("q")
    system.add_lumped_port("p")

    system.gf_model.add_macro("ell", "1.0")
    system.gf_model.add_macro("osc_mass", "1.5")
    system.gf_model.add_macro("spring_k", "5.0")
    system.gf_model.add_macro("damp_c", "0.25")
    system.gf_model.add_macro("force", "0.2*sin(2*pi*t)")

    for brick in (
        S.Brick("M_q", "q*Test_q/ell", [1], dt=True, position="flow"),
        S.Brick("M_p", "p*Test_p/ell", [1], dt=True, position="flow"),
        S.Brick("J_qp", "e_p*Test_q/ell", [1], position="effort"),
        S.Brick("J_pq", "-e_q*Test_p/ell", [1], position="effort"),
        S.Brick("Damping", "-damp_c*e_p*Test_p/ell", [1], position="effort"),
        S.Brick("External_force", "force*Test_p/ell", [1],
                linear=False, position="effort"),
        S.Brick("Spring_law", "(spring_k*q-e_q)*Test_e_q/ell", [1]),
        S.Brick("Momentum_law", "(p/osc_mass-e_p)*Test_e_p/ell", [1]),
    ):
        system.add_brick(brick)

    system.set_initial_value("q", 0.4)
    system.set_initial_value("p", 0.0)
    system.set_time_scheme(
        ts_type="bdf", ts_bdf_order=2,
        ts_exact_final_time="matchstep", ts_adapt_type="none",
        ksp_type="gmres", pc_type="lu",
        pc_factor_mat_solver_type="superlu",
        t_0=0.0, t_f=0.10, dt=0.01, dt_save=0.01,
        init_step=True,
    )
    system.solve()

    # SCRIMP writes its own files under ./outputs/ (workspace cwd).
    # Also export a simple CSV so new users can immediately inspect data.
    times = np.asarray(system.solution["t"], dtype=float).ravel()
    displacement = np.asarray(system.get_solution("q"), dtype=float).ravel()
    momentum = np.asarray(system.get_solution("p"), dtype=float).ravel()
    if not (len(times) == len(displacement) == len(momentum)):
        raise RuntimeError("Unexpected time-series dimensions")
    # Energy is tracked explicitly: this *driven, damped* oscillator
    # exchanges power with the exterior and dissipates energy, so H
    # is not expected to remain constant. Hamiltonian = k*q^2/2+p^2/(2m).
    energy = 0.5 * 5.0 * displacement**2 + momentum**2 / (2 * 1.5)
    output = Path("outputs/oscillator.csv")
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(output, np.column_stack((times, displacement, momentum, energy)),
               delimiter=",", header="time,q,p,H", comments="")
    print(f"Saved {len(times)} time points (including H) to {output.resolve()}")


if __name__ == "__main__":
    main()
