# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Second independent linear DAE qualification: dissipative wave with port.

Construct the *original* examples/wave_dissipative.py equations without
editing the example: replace only its initial momentum, left drive, mesh
resolution and time-scheme for this tiny diagnostic. Audit the mechanical
Hamiltonian, signed boundary work, AND nonnegative damping work.
No experimental init_step monkeypatch: use the production option.
"""

import argparse
import importlib.util
import json
from pathlib import Path
from unittest.mock import patch

from mpi4py import MPI
import numpy as np
import getfem as gf
import scrimp as S
from scrimp.utils.dae_initialization import InconsistentDAEInitialData


def _scalar(value, mpi):
    arr = np.asarray(value).real.ravel()
    if len(arr) != 1:
        raise RuntimeError("Expected scalar energy/power")
    # Native GetFEM asm already globally reduces distributed mesh
    # results; do not MPI-SUM an already global integral.
    maximum = mpi.allreduce(float(arr[0]), op=MPI.MAX)
    minimum = mpi.allreduce(float(arr[0]), op=MPI.MIN)
    if abs(maximum - minimum) > 1e-10 * max(1., abs(maximum)):
        raise RuntimeError(
            "GetFEM global integral differs across MPI ranks"
        )
    return float(maximum)


def _build(h, dt, steps, *, use_default=False):
    # Import the actual example, without copying/re-implementing the
    # viscous constitutive port or the underlying pH weak bricks.
    import_path = (
        Path(__file__).resolve().parents[1]
        / "examples/wave_dissipative.py"
    )
    spec = importlib.util.spec_from_file_location(
        "v130_dissipative_original", import_path
    )
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    original_domain = S.Domain
    original_set_initial = S.DPHS.set_initial_value
    original_set_control = S.DPHS.set_control
    original_set_time = S.DPHS.set_time_scheme
    captured = []

    def domain(name, settings):
        if name == "Rectangle":
            settings = {**settings, "h": h}
        return original_domain(name, settings)

    def initial(self, name, expression):
        if name == "p":
            captured.append(("p",expression))
            expression = "x*(3-x)"
        return original_set_initial(self, name, expression)

    def control(self, name, expression):
        if name == "Boundary control (left)":
            expression = "0.08*sin(20*t)"
        return original_set_control(self, name, expression)

    def scheme(self, **kwargs):
        kwargs.update(
            ts_type="cn", ts_adapt_type="none",
            ts_exact_final_time="matchstep",
            t_f=dt*steps, dt=dt, dt_save=dt,
            init_step=True,
            init_algebraic_atol=1.e-10,
            init_algebraic_rtol=1.e-8,
            init_algebraic_accuracy=1.e-11,
            init_algebraic_max_passes=3,
            ksp_type="preonly", pc_type="lu",
            pc_factor_mat_solver_type="mumps",
        )
        if use_default:
            # The published example is intentionally pinned to
            # historical pseudo for compatibility. Remove that
            # inherited choice only in this new-default diagnostic
            # before invoking the production set_time_scheme API.
            inherited = kwargs.pop("init_step_ts_type", None)
            if inherited not in (None, "pseudo"):
                raise RuntimeError(
                    "Unexpected inherited dissipative-wave strategy"
                )
        else:
            kwargs["init_step_ts_type"] = "algebraic"
        # Distinguish caller kwargs from SCRIMP's materialized PETSc
        # default. PETSc.Options.hasName() is *not* an explicitness
        # signal: it reports both explicit and inserted defaults.
        self._v130_default_probe_explicit_keyword = (
            "init_step_ts_type" in kwargs
        )
        return original_set_time(self, **kwargs)

    def skip_solve(self):
        return None

    def skip_plot(self, **kwargs):
        return None

    with (
        patch.object(S, "Domain", side_effect=domain),
        patch.object(S.DPHS, "set_initial_value", initial),
        patch.object(S.DPHS, "set_control", control),
        patch.object(S.DPHS, "set_time_scheme", scheme),
        patch.object(S.DPHS, "solve", skip_solve),
        patch.object(S.DPHS, "plot_Hamiltonian", skip_plot),
    ):
        system = original.wave_eq()

    if len(captured) != 1:
        raise RuntimeError(
            "Did not intercept exactly one original dissipative-wave p0"
        )
    if not any(
        p.get_name() == "Damping" for p in system.ports.values()
    ):
        raise RuntimeError("Independent dissipative port was lost")
    explicit_kwarg = getattr(
        system, "_v130_default_probe_explicit_keyword", None
    )
    if explicit_kwarg is None:
        raise RuntimeError("Dissipative test did not intercept time scheme")
    if bool(explicit_kwarg) == bool(use_default):
        raise RuntimeError(
            "Dissipative probe initializer keyword presence contradicts "
            "requested implicit/explicit mode"
        )
    if system.time_scheme.getString(
        "init_step_ts_type", ""
    ).strip().lower() != "algebraic":
        raise RuntimeError(
            "Production PETSc.Options did not select algebraic initializer"
        )
    if not system.time_scheme.hasName("init_step_ts_type"):
        raise RuntimeError(
            "SCRIMP did not materialize the effective default"
        )
    return system


def _getfem_region(port):
    """Return the scalar region expected by GetFEM generic assembly.

    SCRIMP Port(region=None) means the entire interior mesh.
    This mirrors scrimp.port.Port.compute(), which maps None to -1.
    Unlike a boundary control, an interior dissipative port has
    no named region; passing None directly to GetFEM is invalid.
    """
    region = port.get_region()
    if region is None:
        return -1
    if isinstance(region, (int, np.integer)):
        return int(region)
    raise TypeError(
        f"Port {port.get_name()} has a noninteger region "
        f"{region!r}; cannot pass it to gf.asm"
    )


def _physical_audit(system, mpi):
    t = np.asarray(system.solution["t"], dtype=float)
    z = np.asarray(system.solution["z"], dtype=float)
    if t.ndim != 1 or t.size < 3 or z.shape[0] != t.size:
        raise RuntimeError("Incomplete dissipative wave snapshots")
    if not np.all(np.diff(t) > 0):
        raise RuntimeError("Non-monotone main CN time grid")
    model = system.gf_model
    H = np.array([
        _scalar(x, mpi) for x in system.get_Hamiltonian()
    ])
    if len(H) != len(t):
        raise RuntimeError("Dissipative Hamiltonian wrong sample count")
    boundary = list(system.controls.values())
    dissipative = [
        p for p in system.ports.values()
        if p.get_name() == "Damping"
    ]
    if len(boundary) != 4 or len(dissipative) != 1:
        raise RuntimeError("Unexpected dissipative-wave port layout")
    ports = boundary + dissipative
    powers = {p.get_name(): [] for p in ports}
    cumulative_boundary, cumulative_dissipation = [0.], [0.]
    balance, defects = [0.], []
    for k in range(1, len(t)):
        zm = 0.5*(z[k] + z[k-1])
        tm = 0.5*(t[k] + t[k-1])
        model.set_time(float(tm))
        model.to_variables(zm)
        values = {}
        for port in ports:
            exp = port.get_flow() + "*" + port.get_effort()
            # SCRIMP Port(region=None) denotes a volume port covering
            # the entire mesh; Port.compute() maps it to GetFEM
            # region=-1.  Passing None to gf.asm is NOT accepted by
            # GetFEM's integer region argument (the native run failed
            # with "argument 5 ... is not a scalar value").
            # Boundary control ports keep their explicit integer region.
            region = _getfem_region(port)
            val = gf.asm(
                "generic", system.domain._int_method[
                    port.get_mesh_id()
                ], 0, exp, region, model
            )
            values[port.get_name()] = _scalar(val, mpi)
            powers[port.get_name()].append(values[port.get_name()])
        dissipated = values["Damping"]
        if dissipated < -1.e-9:
            raise RuntimeError(
                "Damping port dissipates negative power: "
                f"{dissipated:.9g}"
            )
        dt = float(t[k]-t[k-1])
        cumulative_boundary.append(
            cumulative_boundary[-1] +
            dt*sum(values[p.get_name()] for p in boundary)
        )
        cumulative_dissipation.append(
            cumulative_dissipation[-1] + dt*dissipated
        )
        defects.append(float(
            H[k]-H[k-1] +
            cumulative_boundary[-1]-cumulative_boundary[-2] +
            cumulative_dissipation[-1]-cumulative_dissipation[-2]
        ))
        balance.append(float(
            H[k]-H[0]+cumulative_boundary[-1]+cumulative_dissipation[-1]
        ))
    return {
        "times": t.tolist(),
        "hamiltonian": H.tolist(),
        "midpoint_port_powers": powers,
        "integrated_signed_boundary_work": cumulative_boundary[-1],
        "integrated_dissipation": cumulative_dissipation[-1],
        "absolute_dissipation_work": float(sum(
            abs(x) for x in np.diff(cumulative_dissipation)
        )),
        "absolute_signed_boundary_work": float(sum(
            abs(x) for x in np.diff(cumulative_boundary)
        )),
        "balance_defects": balance,
        "step_defects": defects,
        "max_absolute_balance_defect": float(
            max(abs(x) for x in balance)
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=float, default=.2)
    parser.add_argument("--dt", type=float, default=.002)
    parser.add_argument("--steps", type=int, default=10)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--use-default", action="store_true",
        help="Do not set init_step_ts_type, validate library default",
    )
    args = parser.parse_args()
    if args.mesh <= 0 or args.dt <= 0 or args.steps < 2:
        parser.error("Positive mesh, dt and at least two steps required")
    mpi = MPI.COMM_WORLD
    sys = _build(
        args.mesh, args.dt, args.steps,
        use_default=args.use_default,
    )
    prescribed = np.asarray(
        sys.gf_model.from_variables(), dtype=float
    ).copy()
    states = {
        name: tuple(map(int, sys.gf_model.interval_of_variable(name)))
        for name in sys.states
    }
    sys.solve()
    info = sys.initialization_report
    if not info or info["status"] != "PASS":
        raise RuntimeError("Production dissipative algebraic init failed")
    if not info["accuracy_target_met"]:
        raise RuntimeError(
            "Dissipative initial algebraic residual not refined to 1e-11"
        )
    if sys.last_ts_converged_reason <= 0 or (
        sys.last_ts_steps != args.steps or sys.last_ts_rejections
        or sys.last_ts_snes_failures
    ):
        raise RuntimeError("Dissipative physical CN not fully converged")
    t = np.asarray(sys.solution["t"], dtype=float)
    z = np.asarray(sys.solution["z"], dtype=float)
    wanted = np.arange(args.steps+1)*args.dt
    if not np.allclose(t, wanted, atol=1.e-12, rtol=0):
        raise RuntimeError("Unexpected dissipative CN sample times")
    dynamics = {}
    for name, (start, width) in states.items():
        drift = float(np.linalg.norm(
            z[0,start:start+width] -
            prescribed[start:start+width]
        ))
        dynamics[name] = drift
        if drift != 0:
            raise RuntimeError(
                f"Dissipative algebraic initializer altered {name}"
            )
    physical = _physical_audit(sys, mpi)
    # Report actual defect without presuming the accounting on the
    # independent dissipative port is qualified in advance. Final
    # comparator will apply a strict acceptance gate after MPI parity.
    scale = max(1., max(abs(x) for x in physical["hamiltonian"]))
    tolerance = 1.e-9 + 1.e-5*scale
    report = {
        "schema": "scrimp-v130-dissipative-second-linear-dae-v1",
        "rank_count": mpi.Get_size(),
        "implicit_default_requested": bool(args.use_default),
        "explicit_init_step_ts_type_present": bool(
            sys._v130_default_probe_explicit_keyword
        ),
        "effective_PETSc_initial_strategy": sys.time_scheme.getString(
            "init_step_ts_type", ""
        ),
        "mesh": args.mesh, "dt": args.dt, "steps": args.steps,
        "model": "original examples/wave_dissipative.py weak bricks",
        "initial_p0": "x*(3-x)",
        "left_drive": "0.08*sin(20*t)",
        "initial_dynamic_drift": dynamics,
        "dae_projection": info,
        "main_converged_reason": int(sys.last_ts_converged_reason),
        "physical": physical,
        "power_tolerance": tolerance,
    }
    if mpi.Get_rank() == 0:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        stem = f"dissipative_core_{mpi.Get_size()}r"
        (args.output_dir / f"{stem}.json").write_text(
            json.dumps(report, indent=2)+"\n",
            encoding="utf-8",
        )
        np.savez_compressed(
            args.output_dir / f"{stem}.npz",
            times=t, states=z, initial=prescribed,
        )
        print(
            f"SECOND LINEAR DAE MPI={mpi.Get_size()}: "
            f"initial_a={info['native_final_residual_l2']:.9g}, "
            f"drift={dynamics}, dissipative_work="
            f"{physical['integrated_dissipation']:.9g}, "
            f"boundary_work="
            f"{physical['integrated_signed_boundary_work']:.9g}, "
            f"max_power_defect="
            f"{physical['max_absolute_balance_defect']:.9g}, "
            f"tolerance={tolerance:.9g}",
            flush=True,
        )
    if (
        physical["max_absolute_balance_defect"] > tolerance
        or physical["integrated_dissipation"] < -1e-10
    ):
        raise RuntimeError(
            "Independent dissipative Hamiltonian+port+damping "
            "power identity not qualified; report archived"
        )


if __name__ == "__main__":
    main()
