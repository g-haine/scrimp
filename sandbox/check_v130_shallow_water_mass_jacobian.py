# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Native directional FD audit of the shipped shallow_water.py IJacobian.

Load original example weak forms and initial data, intercept only the
long simulation/output, build its actual GetFEM/PETSc operators and
check F'(z,zd)[v] against J(z,zd,shift)v. No pseudo-TS or main time
integration. This is an independent *nonlinear mass* regression, NOT
nonlinear consistent-DAE initialization (which remains unsupported).

The production default must work with SCRIMP_NONLINEAR_MASS_JACOBIAN
UNSET; explicit legacy is a negative control on the same test state.
"""

import argparse
import importlib.util
import json
import os
from pathlib import Path
from unittest.mock import patch

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import scrimp as S


def _build_original(mesh):
    path = Path(__file__).resolve().parents[1] / "examples/shallow_water.py"
    spec = importlib.util.spec_from_file_location(
        "scrimp_v130_original_shallow_water", path
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    original_domain = S.Domain

    def domain(name, options):
        if name == "Rectangle":
            options = {**options, "h": mesh}
        return original_domain(name, options)

    def no_run(self):
        return None

    def no_plot(self, **kwargs):
        return None

    def no_export(self, *args, **kwargs):
        return None

    # Keep every published brick, constitutive expression, initial
    # state and boundary control unchanged. Only omit the 5000-step
    # physical solve and output from the diagnostic setup.
    with (
        patch.object(S, "Domain", side_effect=domain),
        patch.object(S.DPHS, "solve", no_run),
        patch.object(S.DPHS, "plot_Hamiltonian", no_plot),
        patch.object(S.DPHS, "export_to_pv", no_export),
    ):
        system = module.shallow_water_eq()
    if not system.linear_mass and not system.linear_stiffness:
        # set_linear_flags is called only after builder returns.
        pass
    system.set_linear_flags()
    if system.linear_mass:
        raise RuntimeError(
            "Original shipped shallow water is no longer state-mass-dependent"
        )
    # Reproduce the operator setup without creating any TS time step.
    system.allocate_memory()
    system.disable_all_bricks()
    system.gf_model.set_time(0.0)
    system.gf_model.first_iter()
    system.assemble_mass()
    system.assemble_nl_mass()
    system.assemble_stiffness()
    system.assemble_nl_stiffness()
    system.assemble_rhs()
    system.J = system.tangent_stiffness.duplicate()
    system.J.assemble()
    system.F = system.rhs.duplicate()
    system.F.assemble()
    system.buffer = system.rhs.duplicate()
    system.buffer.assemble()
    return system


class _AuditSNES:
    def __init__(self, ksp):
        self.ksp = ksp

    def getKSP(self):
        return self.ksp


class _AuditTS:
    """Expose only the TS interface needed by the real IJacobian.

    No time integration, SNES nonlinear solve or MUMPS factorization.
    The auxiliary KSP uses PC.NONE; it exists only so IJacobian
    can inspect the PETSc PC type and exercise the production path.
    """
    def __init__(self, ksp):
        self.snes = _AuditSNES(ksp)

    def getSNES(self):
        return self.snes

    def getStepNumber(self):
        return 0


def _distributed_vec(template, global_data):
    v = template.duplicate()
    lo, hi = v.getOwnershipRange()
    v.getArray()[:] = global_data[lo:hi]
    return v


def _directional_check(system, ts, z, zd, direction, shift, epsilon):
    m = system.J
    system.IJacobian(ts, 0., z, zd, shift, m, m)
    jv = z.duplicate()
    zplus, zminus = z.copy(), z.copy()
    zdplus, zdminus = zd.copy(), zd.copy()
    fplus, fminus = z.duplicate(), z.duplicate()
    diff = z.duplicate()
    try:
        m.mult(direction, jv)
        zplus.axpy(epsilon, direction)
        zminus.axpy(-epsilon, direction)
        zdplus.axpy(epsilon*shift, direction)
        zdminus.axpy(-epsilon*shift, direction)
        system.IFunction(None, 0., zplus, zdplus, fplus)
        system.IFunction(None, 0., zminus, zdminus, fminus)
        fplus.axpy(-1., fminus)
        fplus.scale(0.5/epsilon)
        fplus.copy(diff)
        diff.axpy(-1., jv)
        fd_norm, jac_norm = float(fplus.norm()), float(jv.norm())
        return {
            "shift": shift,
            "zd_norm": float(zd.norm()),
            "finite_difference_norm": fd_norm,
            "jacobian_action_norm": jac_norm,
            "relative_error": float(diff.norm()) / max(
                fd_norm, jac_norm, 1e-14
            ),
        }
    finally:
        for obj in (
            jv,zplus,zminus,zdplus,zdminus,fplus,fminus,diff
        ):
            obj.destroy()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mesh", type=float, default=.2)
    ap.add_argument("--epsilon", type=float, default=1e-6)
    ap.add_argument("--state-perturb", type=float, default=.01)
    ap.add_argument(
        "--mode", choices=("default", "legacy"), default="default"
    )
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    if args.mesh <= 0 or args.epsilon <= 0 or args.state_perturb < 0:
        ap.error("positive mesh/epsilon and nonnegative state perturb")
    env_mode = os.environ.get("SCRIMP_NONLINEAR_MASS_JACOBIAN")
    if args.mode == "default" and env_mode is not None:
        raise RuntimeError(
            "A default-Jacobian test must leave "
            "SCRIMP_NONLINEAR_MASS_JACOBIAN unset"
        )
    if args.mode == "legacy" and env_mode != "legacy":
        raise RuntimeError(
            "Historical negative control must explicitly select legacy"
        )
    mpi = MPI.COMM_WORLD
    system = _build_original(args.mesh)
    n = int(system.gf_model.nbdof())
    origin = np.asarray(system.gf_model.from_variables()).real.ravel()
    if origin.shape != (n,):
        raise RuntimeError("Wrong original model initial state length")
    rng = np.random.default_rng(513)
    vdata = rng.normal(size=n)
    vdata /= np.linalg.norm(vdata)
    wdata = rng.normal(size=n)
    wdata /= np.linalg.norm(wdata)
    zdata = origin + args.state_perturb * max(
        1., np.linalg.norm(origin)
    ) * vdata
    zddata = .05 * max(1., np.linalg.norm(origin)) * wdata
    template = system.J.createVecRight()
    z = _distributed_vec(template, zdata)
    zd = _distributed_vec(template, zddata)
    zero = template.duplicate()
    zero.set(0.)
    v = _distributed_vec(template, vdata)
    ksp = PETSc.KSP().create(comm=PETSc.COMM_WORLD)
    ksp.getPC().setType(PETSc.PC.Type.NONE)
    audit_ts = _AuditTS(ksp)
    try:
        outcomes = []
        for derivative in (zero, zd):
            for shift in (0., 1.):
                outcomes.append(_directional_check(
                    system, audit_ts, z, derivative, v, shift, args.epsilon
                ))
        residual_max = max(row["relative_error"] for row in outcomes)
        report = {
            "schema":"scrimp-v130-shipped-shallow-water-mass-jacobian-v1",
            "model":"original examples/shallow_water.py weak bricks",
            "mode":args.mode,
            "ranks":mpi.Get_size(),
            "dofs":n,
            "initial_state_perturbation":args.state_perturb,
            "epsilon":args.epsilon,
            "nonlinear_mass":not system.linear_mass,
            "mass_mode_env":env_mode,
            "main_TS_executed":False,
            "independent_time_derivative":True,
            "outcomes":outcomes,
            "maximum_directional_relative_error":residual_max,
        }
        if mpi.Get_rank()==0:
            args.output_dir.mkdir(parents=True,exist_ok=True)
            (args.output_dir/f"shallow_{args.mode}_{mpi.Get_size()}r.json").write_text(
                json.dumps(report, indent=2)+"\n",encoding="utf-8"
            )
            print(
                f"SHALLOW MASS JACOBIAN mode={args.mode} MPI={mpi.Get_size()} "
                f"DOFs={n} max_rel={residual_max:.9g} "
                f"nonzero_zd_rel={[r['relative_error'] for r in outcomes[2:]]}",
                flush=True,
            )
        if not np.isfinite(residual_max):
            raise FloatingPointError("Nonfinite production Jacobian error")
    finally:
        for obj in (v,zero,zd,z,template,ksp):
            obj.destroy()


if __name__=="__main__":
    main()
