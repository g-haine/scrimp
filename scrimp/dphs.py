# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             dphs.py
- authors:          Giuseppe Ferraro, Ghislain Haine, Florian Monteghetti
- date:             22 nov. 2022
- brief:            class for distributed port-Hamiltonian system
"""

from scrimp.utils.linalg import convert_gmm_to_petsc, extract_gmm_to_scipy
from scrimp.utils.profiling import TimingProfiler
from scrimp.utils.initial_state import publish_accepted_pseudo_state
from scrimp.utils.dae_initialization import initialize_linear_algebraic
from scrimp.utils.fe_initialization import project_fe_initial_value
from scrimp.utils.dae_checkpoint import (
    model_provenance as _checkpoint_model_provenance,
    collective_read as _checkpoint_collective_read,
    collective_write as _checkpoint_collective_write,
)
from scrimp.utils.nonlinear_affine_dae import (
    initialize_nonlinear_affine_algebraic,
)
from scrimp.utils.matrix_memory import (
    alias_readonly_linear_tangent, add_mass_hessian_fd_reused,
)
from scrimp.utils.phase_memory_audit import PhaseMemoryAudit
from scrimp.utils.power_integration import integrate_port_power
from scrimp.utils.async_binary_trajectory import (
    RankSender, LazyTrajectory, wait_for_manifest,
)
from scrimp.hamiltonian import Hamiltonian
from scrimp.brick import Brick
from scrimp.control import Control_Port
from scrimp.fem import FEM
from scrimp.port import Parameter, Port
from scrimp.costate import CoState
from scrimp.state import State
from scrimp.domain import Domain
from petsc4py import PETSc
import getfem as gf
import logging
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np
import os
import re
from mpi4py import MPI

import petsc4py
import sys
import time
import gc

petsc4py.init(sys.argv)
comm = PETSc.COMM_WORLD
rank = comm.getRank()
max_rank = comm.getSize()


def _mpi_trace(stage):
    """Emit per-rank milestones without adding MPI collectives/barriers."""
    if os.environ.get("SCRIMP_MPI_TRACE") == "1":
        print(f"[SCRIMP MPI {rank}/{max_rank}] {stage}",
              file=sys.stderr, flush=True)


def _set_main_ts_default_exact_final_time(ts):
    """Default main time integration to the physically integrated t_f.

    PETSc's ts_exact_final_time option can still override this choice
    when a caller explicitly requests legacy STEPOVER semantics.
    """
    ts.setExactFinalTime(PETSc.TS.ExactFinalTime.MATCHSTEP)


def _select_saved_time_indices_through_final_time(times, t_final):
    """Select time metadata and archive indices without stacking states."""
    t = np.asarray(times, dtype=float)
    if not np.isfinite(t_final):
        raise ValueError("Final time must be finite")
    ulps = 32.0 * abs(np.spacing(t_final if t_final else 1.0))
    order = np.argsort(t)
    kept = order[t[order] <= t_final + ulps]
    selected = t[kept].copy()
    selected[np.abs(selected - t_final) <= ulps] = t_final
    return selected, kept


def _select_saved_solution_through_final_time(times, states, t_final):
    """Keep accepted solution snapshots up to t_f with ULP tolerance.

    Floating-point step accumulation may give t_final+1 ULP when the
    mathematical endpoint is exactly t_final. A strict <= t_final
    mask drops that correct final snapshot (MCSS2024 reproduced:
    PETSc reached 0.00030000000000000003, t_f stored as 0.0003).
    Tolerance is 32 representable spacings *at t_final*, rather
    than an absolute 1e-12 that could be too coarse for tiny horizons.
    Truly STEPOVER steps beyond the requested horizon remain excluded.
    """
    # Filter time indices before materializing any state array.
    t, kept = _select_saved_time_indices_through_final_time(times, t_final)
    if isinstance(states, np.ndarray):
        z = states[kept]  # One indexed copy, same output dtype/shape.
    elif len(kept):
        # Exactly one stacked allocation; do not stack rejected states.
        z = np.stack([states[int(i)] for i in kept], axis=0)
    elif len(states):
        # Match historical (0, ...) shape even if every state is rejected.
        sample = np.asarray(states[0])
        z = np.empty((0,) + sample.shape, dtype=sample.dtype)
    else:
        z = np.asarray(states)
    return t, z


# GetFEM promises a globally summed RHS, replicated on all ranks, NOT
# a separate piece to be summed again in PETSc:
# https://getfem.org/userdoc/parallel.html
#
# A full replicated RHS is normal and does NOT justify division by P.
# On the original local GetFEM 5.4.2 MPI-2 build, however, actual
# integrals were P times the sequential reference at P=1,2,4, while
# the full raw tangent matrices were replicated, rank for rank.
# build_rhs and build_all produce the same discrepancy. The underlying
# partition/reduction or build issue is DEFERRED in GitHub issue #23;
# the exact GetFEM implementation fault has not been established.
# Avoid assuming ALL GetFEM MPI-2 binaries have that behavior. Calibrate
# the installed generic-assembly RHS once against an exactly known
# integral and normalize only if it really is duplicated.
_getfem_rhs_calibration = None


def _getfem_mpi_rhs_scale():
    """Calibrate MPI reduction against the exact integral of u on [0,1].

    A GetFEM reference model with u=1 and weak form u*Test_u must
    return a negative RHS whose entry sum is -1, regardless of P.
    Detect either correct global summation (1) or P-fold duplicated
    summation (P). Do not silently normalize unexpected values.
    This probes *generic nonlinear bricks* only: custom/explicit
    RHS bricks can have different MPI conventions.
    """
    global _getfem_rhs_calibration
    if max_rank == 1 or int(gf.util_mpi_parallelism_level()) != 2:
        return 1.0
    if _getfem_rhs_calibration is not None:
        return _getfem_rhs_calibration

    # All ranks must execute GetFEM's MPI-parallel model assembly.
    mesh = gf.Mesh("cartesian", np.linspace(0.0, 1.0, 17))
    mf = gf.MeshFem(mesh, 1)
    mf.set_fem(gf.Fem("FEM_PK(1,1)"))
    mim = gf.MeshIm(mesh, gf.Integ("IM_GAUSS1D(5)"))
    probe = gf.Model("real")
    probe.add_fem_variable("rhs_probe", mf)
    probe.add_nonlinear_term(mim, "rhs_probe*Test_rhs_probe")
    probe.set_variable("rhs_probe", np.ones(probe.nbdof()))
    probe.assembly(option="build_rhs")
    factor = -float(np.real(np.asarray(probe.rhs()).sum()))
    # Respect the PETSc communicator; verify all processes independently
    # see the same GetFEM reduction factor before accepting it.
    factors = comm.tompi4py().allgather(factor)
    if not np.allclose(factors, factor, rtol=1e-10, atol=1e-10):
        raise RuntimeError(
            "Inconsistent GetFEM RHS normalization across MPI ranks: "
            f"{factors}"
        )
    if np.isclose(factor, 1.0, rtol=1e-8, atol=1e-8):
        scale = 1.0
    elif np.isclose(factor, max_rank, rtol=1e-8, atol=1e-8):
        scale = 1.0 / max_rank
    else:
        raise RuntimeError(
            "Unexpected GetFEM MPI RHS integral for u=1 on [0,1]: "
            f"{factor:.12g} (expected 1 or {max_rank}). "
            "Set SCRIMP_GETFEM_RHS_MODE=raw only after verifying "
            "the installed GetFEM assembly."
        )
    _getfem_rhs_calibration = scale
    _mpi_trace(
        f"GetFEM RHS calibration: measured={factor:.12g}, "
        f"mpi_ranks={max_rank}, scale={scale:.12g}"
    )
    return scale


def _getfem_source_rhs_mode():
    """Opt-in GetFEM MPI source-brick treatment.

    The tested GetFEM MPI-2 build duplicates generic nonlinear RHS
    contributions P times, but add_source_term_brick is already globally
    reduced. Keep the old uniform treatment unless explicitly requested.
    """
    mode = os.environ.get(
        "SCRIMP_GETFEM_SOURCE_RHS_MODE", "match_generic"
    ).strip().lower()
    if mode not in ("match_generic", "raw"):
        raise ValueError(
            "SCRIMP_GETFEM_SOURCE_RHS_MODE must be match_generic or raw"
        )
    return mode


def _get_mumps_reset_policy(all_bricks_linear=False):
    """Select the MUMPS reset policy, honoring explicit user overrides.

    With no override: on_change for ALL models (including nominally
    linear ones). Exact CSR sparsity and coefficient comparisons are
    necessary because linear does not imply time-invariant, and TS shift
    can change under variable-step integration.
    The legacy boolean flag is still accepted for old benchmarks.
    An explicit policy ALWAYS wins over automatic linearity detection.
    """
    policy = os.environ.get("SCRIMP_MUMPS_RESET_POLICY")
    if policy is None:
        legacy = os.environ.get("SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN")
        if legacy is None:
            return "on_change"
        if legacy == "1":
            return "always"
        if legacy == "0":
            return "never"
        raise ValueError(
            "SCRIMP_MUMPS_RESET_PC_ON_JACOBIAN must be 0 or 1"
        )
    policy = policy.strip().lower()
    if policy not in ("on_change", "always", "never"):
        raise ValueError(
            "SCRIMP_MUMPS_RESET_POLICY must be on_change, always or never"
        )
    return policy


# Known GWFL time, geometric symbols and elementary functions only.
# Any unknown identifier is deliberately classified as potentially
# state-dependent: custom macros may hide references to GetFEM variables.
# Do NOT infer that all "source" bricks are independent of the state.
_PRESCRIBED_CONTROL_SYMBOLS = frozenset({
    "t", "x", "y", "z", "X", "Normal", "Tangent", "pi", "Pi",
    "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
    "sinh", "cosh", "tanh", "exp", "log", "sqrt", "pow",
    "min", "max", "abs", "sign", "floor", "ceil",
    "Norm", "Trace", "Grad", "Div", "Det", "Id", "Sym",
})


def _control_expression_is_prescribed(expression, model_variables):
    """Certify a simple open-loop GWFL control by its lexical symbols.

    All model unknowns take precedence over the safe symbol allowlist.
    Unknown functions/macros/data yield False (conservative nonlinear
    feedback path). This deliberately does NOT attempt to parse or
    expand arbitrary GetFEM user-defined macros.
    """
    if not isinstance(expression, str) or not expression.strip():
        return False
    # Scientific notation (e.g. 1e-8, 2.5E+4) must not be mistaken
    # for a model symbol named 'e' or 'E'. Remove numeric literals
    # before extracting GWFL identifiers, without evaluating code.
    without_exponents = re.sub(
        r"(?<![\w.])(?:\d+(?:\.\d*)?|\.\d+)[eE][+-]?\d+\b",
        "0", expression,
    )
    symbols = set(re.findall(r"\b[A-Za-z_][A-Za-z0-9_]*\b", without_exponents))
    # Match complete identifiers, so Y and Y_n cannot be confused.
    unknowns = symbols.difference(_PRESCRIBED_CONTROL_SYMBOLS)
    return not symbols.intersection(model_variables) and not unknowns


def _control_expression_is_direct_linear_feedback(expression, model_variables):
    """Certify only signed *bare unknowns* as constant linear feedback.

    For example `Y_w`, `+Y_T`, and `-Y_T` yield a GetFEM
    bilinear weak term `Y_w*Test_U_T` with a fixed Jacobian.
    Products, arbitrary coefficients, macros, and functions must
    stay on the nonlinear GetFEM route. This fail-closed test
    deliberately does not claim to recognize every affine law.
    """
    if not isinstance(expression, str):
        return False
    matched = re.fullmatch(
        r"\s*[+-]?\s*([A-Za-z_][A-Za-z0-9_]*)\s*", expression
    )
    return bool(matched and matched.group(1) in model_variables)


def _all_bricks_are_linear(bricks):
    """Only skip CSR checks for provably affine, fixed-graph systems.

    A GetFEM source brick has no nonlinear Jacobian contribution,
    but we only trust sources *created by set_control* and certified
    to be state-independent. User-supplied or unknown source bricks
    remain conservative. Bare signed model-variable feedback is
    assembled as an exact linear GetFEM tangent; genuinely nonlinear
    or unknown feedback remains nonlinear with its Newton derivative.
    """
    return bool(bricks) and all(
        brick.get_linear()
        or (
            brick.get_position() == "source"
            and getattr(brick, "_scrimp_prescribed_control_source", False)
        )
        for brick in bricks.values()
    )


def _local_stored_csr_pattern(matrix):
    """Snapshot exact stored row pointers/columns of the owned PETSc rows.

    This is the same comparison as the validated standalone PETSc
    on_change probes. The owned-row interval is also recorded, so a
    repartition can never masquerade as unchanged local CSR.
    """
    rowptr, columns, _values = matrix.getValuesCSR()
    return (
        tuple(int(x) for x in matrix.getOwnershipRange()),
        np.array(rowptr, copy=True),
        np.array(columns, copy=True),
    )


def _local_csr_pattern_changed(before, after):
    return (
        before[0] != after[0]
        or not np.array_equal(before[1], after[1])
        or not np.array_equal(before[2], after[2])
    )


def _local_stored_csr_snapshot(matrix):
    """Exact CSR graph AND numeric values for owned PETSc rows.

    Do not use fingerprints or mathematical linearity assumptions:
    a changed coefficient or TS shift requires a fresh numeric LU.
    The snapshot is per rank and includes ownership to catch repartition.
    """
    rowptr, columns, values = matrix.getValuesCSR()
    return (
        tuple(int(x) for x in matrix.getOwnershipRange()),
        np.array(rowptr, copy=True),
        np.array(columns, copy=True),
        np.array(values, copy=True),
    )


def _local_csr_change_flags(before, after):
    """Bit 0: graph changed; bit 1: numeric values changed.

    A changed graph necessarily implies a changed operator. Exact
    comparisons deliberately treat NaN as changed, never as reusable.
    """
    if before is None:
        return 3
    pattern_changed = _local_csr_pattern_changed(before, after)
    numeric_changed = (
        pattern_changed
        or not np.array_equal(before[3], after[3])
    )
    return int(pattern_changed) | (int(numeric_changed) << 1)


import scrimp.utils.config
outputs_path = scrimp.utils.config.outputs_path

class DPHS:
    """A generic class handling distributed pHs using the GetFEM tools

    This is a wrapper in order to simplify the coding process

    Access to fine tunings is preserved as much as possible
    """

    def __init__(self, basis_field="real"):
        """Constructor of the `distributed port-Hamiltonian system` dphs class.

        Args:
            basis_field (str): basis field for unknowns (must be `real` or `complex`)
        """

        #: The `domain` of a dphs is an object that handle mesh(es) and dict of regions with getfem indices (for each mesh), useful to define `bricks` (i.e. forms) in the getfem syntax
        self.domain = None
        #: Opt-in, rank-local timing collector (SCRIMP_PROFILE=1).
        self.profiler = TimingProfiler()
        self._profile_phase = "setup"
        #: True-residual diagnostic for the optional algebraic DAE
        #: initializer. Reset for each solve; None for legacy modes.
        self.initialization_report = None
        #: Clear and init `time_scheme` member, which embed PETSc TS options database
        self.get_cleared_TS_options()
        #: The dict of `states`, store many infos for display()
        self.states = dict()
        #: The dict of `costates`, store many infos for display()
        self.costates = dict()
        #: The dict of `ports`, store many infos for display()
        self.ports = dict()
        #: Spatially global scalar variables (one GetFEM DOF each)
        self.lumped_variables = set()
        #: The dict of `bricks`, associating getfem `bricks` to petsc matrices obtained by the PFEM, store many infos for display()
        self.bricks = dict()
        #: The dict of `controls`, collecting information about control ports. Also appear in `ports` entries
        self.controls = dict()
        #: The `Hamiltonian` of a dphs is a list of dict containing several useful information for each term
        self.hamiltonian = Hamiltonian("Hamiltonian")
        #: To check if the powers have been computed
        self.powers_computed = False
        #: Tangent (non-linear + linear) mass matrix of the system in PETSc CSR format
        self.tangent_mass = PETSc.Mat().create(comm=comm)
        #: Non-linear mass matrix of the system in PETSc CSR format
        self.nl_mass = PETSc.Mat().create(comm=comm)
        #: Linear mass matrix of the system in PETSc CSR format
        self.mass = PETSc.Mat().create(comm=comm)
        #: Tangent (non-linear + linear) stiffness matrix of the system in PETSc CSR format
        self.tangent_stiffness = PETSc.Mat().create(comm=comm)
        #: Non-linear stiffness matrix of the system in PETSc CSR format
        self.nl_stiffness = PETSc.Mat().create(comm=comm)
        #: Linear stiffness matrix of the system in PETSc CSR format
        self.stiffness = PETSc.Mat().create(comm=comm)
        #: rhs of the system in PETSc Vec
        self.rhs = PETSc.Vec().create(comm=comm)
        #: To check if the initial values have been set before time-resolution
        self.initial_value_set = dict()
        #: To check if the PETSc TS time-integration parameters have been set before time-integration
        self.time_scheme["isset"] = False
        #: For monitoring time in TS resolution
        self.ts_start = 0
        #: Will contain both time t and solution z
        self.solution = dict()
        #: Time t where the solution have been saved
        self.solution["t"] = list()
        #: Solution z at time t
        self.solution["z"] = list()
        #: Only ARCHIVED states change backend; PETSc/TS/Newton state does not.
        self._trajectory_backend = "memory"
        self._trajectory_sender = None
        self._trajectory_directory = None
        #: Keep historical matrix layout until compact mode is validated.
        self._matrix_memory_mode = "legacy"
        self._mass_hessian_fd_probe = None
        self._mass_hessian_fd_matrix = None
        self._mass_hessian_fd_probe_allocations = 0
        self._mass_hessian_fd_matrix_allocations = 0
        self._mass_hessian_fd_matrix_reuses = 0
        # A completely explicit nonlinear brick is assembled only at
        # the previously accepted TS state/time, never at a Newton trial.
        # The PETSc-owned RHS stays fixed during each implicit step.
        self._explicit_rhs = None
        self._explicit_rhs_ready = False
        self._explicit_rhs_time = None
        #: To check if the system has been solved
        self.solve_done = False
        # A verified full accepted state pending a cold restart. It is
        # applied only AFTER GetFEM/PETSc matrix allocation and assembly,
        # BEFORE the first main physical PETSc TS step. No BDF history.
        self._checkpoint_pending = None
        self.last_checkpoint_report = None
        #: To stop TS integration by keeping already computed timesteps if one step fails
        self.stop_TS = False
        #: To speed-up IFunction calls for linear systems
        self.linear_mass = True
        self.linear_stiffness = True
        #: A PETSc Vec for residual computation
        self.F = PETSc.Vec().create(comm=comm)
        #: A PETSc Mat for Jacobian computation
        self.J = PETSc.Mat().create(comm=comm)
        #: A PETSc Vec buffering computation
        self.buffer = PETSc.Mat().create(comm=comm)
        #: A getfem `Model` object that is use as core for the dphs
        self.gf_model = gf.Model(basis_field)
        #: Says to getfem that the `Model` is time-dependent
        self.gf_model.add_initialized_data("t", 0.0, sizes=1)

        # HACK: Macros x, y and z are not available for source term otherwise (why?)
        # Relative problem: source term does not accept numpy functions for the moment (why?)
        self.gf_model.add_macro("x", "X(1)")
        self.gf_model.add_macro("y", "X(2)")
        self.gf_model.add_macro("z", "X(3)")

        if rank == 0:
            logging.info(f"A model with {basis_field} unknowns has been initialized")

    def _checkpoint_provenance_collective(self):
        """Reject a missing or divergent full GetFEM model on all MPI ranks."""
        mpi = comm.tompi4py()
        provenance = None
        local_error = None
        try:
            provenance = _checkpoint_model_provenance(self)
        except Exception as exc:
            local_error = f"{type(exc).__name__}: {exc}"
        errors = mpi.allgather(local_error)
        if any(error is not None for error in errors):
            raise ValueError(
                "Complete native GetFEM checkpoint provenance cannot be "
                "established collectively: " + repr(errors)
            )
        hashes = mpi.allgather(provenance["fingerprint_sha256"])
        if any(value != hashes[0] for value in hashes):
            raise ValueError(
                "MPI ranks disagree on physical mesh/FEM/DAE model fingerprint"
            )
        return provenance

    def save_checkpoint(self, path, *, extra=None):
        """Save a complete accepted final state for a later COLD restart.

        Only supported for completed, accepted PETSc TSCN simulations.
        Includes differential and algebraic GetFEM variables and the
        absolute physical time. The file is atomic and never overwritten.
        This is NOT a snapshot of PETSc TS/KSP/MUMPS or BDF history.
        This operation is collective: call it on every MPI rank.
        """
        mpi = comm.tompi4py()
        error = None
        state = None
        final_time = None
        try:
            if not self.solve_done:
                raise ValueError("Cannot checkpoint before a completed solve")
            if int(self.last_ts_converged_reason) <= 0:
                raise ValueError("Cannot checkpoint a nonconverged PETSc TS")
            if int(self.last_ts_rejections) or int(self.last_ts_snes_failures):
                raise ValueError(
                    "Cold checkpoint MVP refuses TS with rejected/failed steps"
                )
            if self.time_scheme.getString("ts_type", "").lower() != "cn":
                raise ValueError("Only PETSc TSCN cold state restart is supported")
            final_time = float(self.last_ts_time)
            horizon = float(self.time_scheme["t_f"])
            if not np.isfinite(final_time) or abs(final_time-horizon) > (
                1e-12*max(1.,abs(horizon))
            ):
                raise ValueError("Actual accepted TS time differs from t_f")
            if not len(self.solution["t"]) or not len(self.solution["z"]):
                raise ValueError("No complete accepted final state was archived")
            if abs(float(self.solution["t"][-1])-final_time) > (
                1e-12*max(1.,abs(final_time))
            ):
                raise ValueError(
                    "Archived final state is not the accepted PETSc TS endpoint"
                )
            state = np.asarray(self.solution["z"][-1], dtype=float).copy()
            if state.ndim != 1 or state.size != int(self.gf_model.nbdof()):
                raise ValueError("Archived final state lacks full FE coefficients")
            if not np.all(np.isfinite(state)):
                raise ValueError("Archived final FE state is not finite")
        except Exception as exc:
            error = f"{type(exc).__name__}: {exc}"
        errors = mpi.allgather(error)
        if any(item is not None for item in errors):
            raise ValueError("Collective accepted-state checkpoint rejected: "+
                             repr(errors))
        provenance = self._checkpoint_provenance_collective()
        return _checkpoint_collective_write(
            path, mpi=mpi, z=state, t=final_time, provenance=provenance,
            ts_type="cn", accepted=True, extra=extra,
        )

    def load_checkpoint(self, path):
        """Stage a verified full DAE state for a new, cold PETSc CN solve.

        Call after building the identical FE model and before solve().
        Sets absolute physical t_0 from checkpoint, forces init_step=false
        (no pseudo/BDF or algebraic consistency reinitialization), then
        rechecks original full GetFEM algebraic equations after assembly.
        The future physical boundary law must be the SAME as the saved
        model. Collective on all MPI ranks. Returns checkpoint metadata.
        """
        if self.solve_done or self._checkpoint_pending is not None:
            raise ValueError(
                "Loading an accepted-state checkpoint requires a fresh DPHS"
            )
        if not self.time_scheme.getBool("isset"):
            raise ValueError("Set a physical PETSc time scheme before restart")
        if self.time_scheme.getString("ts_type", "").lower() != "cn":
            raise ValueError("Only PETSc TSCN cold state restart is supported")
        provenance = self._checkpoint_provenance_collective()
        mpi = comm.tompi4py()
        t_saved, state, header = _checkpoint_collective_read(
            path, mpi=mpi, current_provenance=provenance, ts_type="cn",
        )
        t_final = float(self.time_scheme["t_f"])
        t_initial = float(self.time_scheme["t_0"])
        if not (np.isfinite(t_final) and t_final > t_saved):
            raise ValueError("Checkpoint must precede future physical t_f")
        if abs(t_initial) > 1e-12 and abs(t_initial-t_saved) > 1e-12:
            raise ValueError(
                "Configured absolute t_0 conflicts with checkpoint time"
            )
        self.set_time_scheme(t_0=t_saved, init_step=False)
        self._checkpoint_pending = (t_saved, state.copy(), header)
        return header

    def _publish_pending_checkpoint_after_assembly(self):
        """Verify native DAE algebraic rows before mutating the initial FE state.

        Called after original physical GetFEM/PETSc assembly but before
        PETSc TS is created. Validates F_a(t*,z*,0), NOT full F(t*,z*,0):
        a moving physical state may have nonzero differential derivatives.
        """
        pending = self._checkpoint_pending
        if pending is None:
            return
        time_saved, state, header = pending
        mpi = comm.tompi4py()
        n = int(self.gf_model.nbdof())
        if len(state) != n or abs(
            float(self.time_scheme["t_0"])-time_saved
        ) > 1e-12:
            raise ValueError(
                "Checkpoint full FE vector or absolute physical clock changed"
            )
        if self.time_scheme.getString("init_step", "true").lower() != "false":
            raise ValueError(
                "Cold checkpoint may not run a pseudo/BDF initializer"
            )
        dynamic = []
        for dynamic_state in self.states.values():
            name = dynamic_state.get_name()
            lo, count = map(int, self.gf_model.interval_of_variable(name))
            dynamic.extend(range(lo, lo + count))
        if not dynamic or len(set(dynamic)) != len(dynamic):
            raise ValueError(
                "Full DAE state has no valid original differential variables"
            )
        algebraic = np.ones(n, dtype=bool)
        algebraic[np.asarray(dynamic, dtype=int)] = False
        vec = self.tangent_stiffness.createVecRight()
        zd = vec.duplicate()
        residual = self.F.duplicate()
        previous_time = float(self.time_scheme["t_0"])
        previous = np.asarray(
            self.gf_model.from_variables(), dtype=float
        ).copy()
        errors = None
        fa_norm = None
        try:
            lo, hi = vec.getOwnershipRange()
            if hi > lo:
                vec.setValues(
                    range(lo, hi), state[lo:hi],
                    addv=PETSc.InsertMode.INSERT_VALUES,
                )
            vec.assemble()
            zd.set(0.0)
            residual.set(0.0)
            self.gf_model.set_time(time_saved)
            self.IFunction(None, time_saved, vec, zd, residual)
            owned = np.asarray(residual.getArray(readonly=True), dtype=float)
            relevant = owned[algebraic[lo:hi]]
            local_square = float(np.dot(relevant, relevant))
            fa_norm = float(np.sqrt(mpi.allreduce(local_square, op=MPI.SUM)))
            if not np.isfinite(fa_norm) or fa_norm >= 1e-9:
                raise ValueError(
                    f"Original GetFEM algebraic checkpoint incompatible: "
                    f"||F_a||={fa_norm:.9g}"
                )
            # No adjustment to the prescribed differential variables:
            # the original complete native state is published bitwise.
            self.gf_model.set_time(time_saved)
            self.gf_model.to_variables(state.copy())
            self.GlobalVec.getArray()[:] = state
            exact = np.array_equal(
                np.asarray(self.gf_model.from_variables(), dtype=float),
                state
            ) and np.array_equal(
                np.asarray(self.GlobalVec.getArray(readonly=True),
                           dtype=float), state
            )
            if not mpi.allreduce(exact, op=MPI.LAND):
                raise ValueError("Cold FE state publication changed coefficients")
        except Exception:
            self.gf_model.to_variables(previous)
            self.gf_model.set_time(previous_time)
            raise
        finally:
            residual.destroy()
            zd.destroy()
            vec.destroy()
        self.last_checkpoint_report = {
            "status": "COMPLETE_DAE_STATE_COLD_RESTORED_IN_DPHS",
            "physical_time": time_saved,
            "original_algebraic_residual_l2": fa_norm,
            "model_sha256": header["model"]["fingerprint_sha256"],
            "no_pseudo_BDF_or_projection": True,
            "all_GetFEM_and_GlobalVec_coefficients_restored_bitwise": True,
            "PETSc_TS_history_restored": False,
        }
        self._checkpoint_pending = None

    def get_profile(self):
        """Collect per-rank timings, call counts and MPI min/mean/max.

        **Collective** when SCRIMP_PROFILE=1: every MPI rank must call
        this method after solve(). Returns None when profiling is off.
        Events are inclusive, including nested GetFEM and PETSc work.
        """
        return self.profiler.report(comm.tompi4py())

    def _profile_ifunction(self, *args):
        with self.profiler.measure(
            f"{self._profile_phase}.petsc.callback.ifunction"
        ):
            return self.IFunction(*args)

    def _profile_ijacobian(self, *args):
        with self.profiler.measure(
            f"{self._profile_phase}.petsc.callback.ijacobian"
        ):
            return self.IJacobian(*args)

    def set_domain(self, domain: Domain):
        """This function sets a domain for the dphs.

        TODO: If not built_in, given from a script 'name.py' or a .geo file with args in the dict 'parameters' should be able to handle several meshes e.g. for interconnections, hence the list type

        Args:
            name (str): id of the domain, either for built in, or user-defined auxiliary script
            parameters (dict): parameters for the construction, either for built in, or user-defined auxiliary script
        """

        self.domain = domain
        if rank == 0:
            logging.info(f"domain: {domain.get_name()} has been set")

    def add_state(self, state: State):
        """This functions adds a state to state dict of the dphs.

        Args:
            state (State): the state
        """

        self.states[state.get_name()] = state

        if rank == 0:
            logging.info(f"state: {state.get_name()} has been added")

    def add_costate(self, costate: CoState):
        """This function adds a costate to the costate dict of the dphs, then defines and adds the dynamical port gathering the couple (State, CoState).

        Args:
            costate (CoState): the costate
        """

        # Set the `costate` in its `state`
        state = costate.get_state()
        state.set_costate(costate)

        # Add the `costate` to the list
        self.costates[costate.get_name()] = costate

        # Define the `port` gathering the `state` and the `costate`
        port = Port(
            state.get_name(),
            state.get_name(),
            costate.get_name(),
            costate.get_kind(),
            state.get_mesh_id(),
            algebraic=False,
            dissipative=False,
            substituted=costate.get_substituted(),
            region=state.get_region(),
        )
        # Add the `port` the the dphs
        self.add_port(port)

        # Set the `port` in the `state` and the `costate`
        state.set_port(port)
        costate.set_port(port)

        if rank == 0:
            logging.info(
                f"costate: {costate.get_name()} has been added to state: {costate.get_state().get_name()}"
            )
            logging.info(
                f"state: {state.get_name()} has new costate: {state.get_costate().get_name()}"
            )

    def add_port(self, port: Port):
        """This function adds a port to the port dict of the dphs

        Args:
            port (Port): the port
        """

        self.ports[port.get_name()] = port
        where = ""
        if port.get_region() is not None:
            where = f"on region {port.get_region()}"

        if rank == 0:
            logging.info(f"port: {port.get_name()} has been added {where}")

    def audit_interface_ownership(self, interface, subdomains, mesh_id=0):
        """Audit an interface against this DPHS's declared ports and bricks.

        The public GetFEM interface region is not modified. The domain
        verifies real two-sided geometry; the ownership audit then uses
        *only* explicit bulk-variable anchors in actual registered weak
        forms. It never silently reorients FEMs or bricks and must not
        be interpreted as a physical/model/DAE qualification.

        Call after defining the ports and bricks (including controls).
        Raises InterfaceOwnershipError for missing/conflicting owners,
        and InterfaceTopologyError for unsupported interface geometry.
        """
        from scrimp.interface_ownership import (
            InterfaceOwnershipError, infer_interface_port_owners,
        )

        if self.domain is None or not self.domain.get_isSet():
            raise ValueError("Domain must be set before interface audit")
        if not isinstance(mesh_id, int) or not 0 <= mesh_id < len(self.domain.get_mesh()):
            raise ValueError(f"Invalid mesh_id for interface: {mesh_id}")
        subdomains = tuple(int(r) for r in subdomains)
        if len(subdomains) != 2 or len(set(subdomains)) != 2:
            raise InterfaceOwnershipError("Two distinct subdomains required")
        interface = int(interface)

        # Prove topology first, including opposing parent-face normals.
        sides = self.domain.resolve_interface(
            interface, subdomains, mesh_id=mesh_id,
        )

        bulk_variables = {}
        interface_variables = {}
        for port in self.ports.values():
            if port.get_mesh_id() != mesh_id:
                continue
            region = port.get_region()
            names = (port.get_flow(), port.get_effort())
            if region in subdomains:
                for name in names:
                    if not name:
                        continue
                    previous = bulk_variables.setdefault(name, int(region))
                    if previous != int(region):
                        raise InterfaceOwnershipError(
                            f"Bulk variable {name} belongs to multiple regions"
                        )
            elif region == interface:
                interface_variables[port.get_name()] = names

        if not interface_variables:
            raise InterfaceOwnershipError(
                f"No declared interface ports on region {interface}"
            )
        forms = (
            (brick.get_regions(), brick.get_form())
            for brick in self.bricks.values()
            if brick.get_mesh_id() == mesh_id
        )
        owners = infer_interface_port_owners(
            interface, interface_variables, bulk_variables, forms,
        )
        if not set(owners.values()).issubset(set(subdomains)):
            raise InterfaceOwnershipError(
                "Interface port owner is not an adjacent subdomain"
            )
        return {
            "interface": interface,
            "mesh_id": mesh_id,
            "subdomains": subdomains,
            "facets_per_side": len(sides.facet_keys),
            "port_owners": owners,
            "topology_verified": True,
            "assembled_with_oriented_faces": False,
        }

    def add_FEM(self, fem: FEM):
        """This function adds a FEM (Finite Element Method) for the variables associated to a port of the dphs

        TODO: handle `tensor-field`

        Args:
            fem (FEM): the FEM to use
        """

        # Check if the `domain` is set
        try:
            assert self.domain.get_isSet()
        except AssertionError as err:
            logging.error("Domain must be set before adding FEM")
            raise err

        # Check if the `port` is in the dict of the dphs
        name_port = fem.get_name()
        try:
            assert name_port in self.ports.keys()
        except AssertionError as err:
            logging.error(f"Port {name_port} does not exist")
            raise err

        # Select the dimension, according to the `kind` of the `port`
        port = self.ports[name_port]
        if port.get_kind() == "scalar-field":
            dim = 1
        elif port.get_kind() == "vector-field":
            dim = self.domain.get_dim()[port.get_mesh_id()]
        else:
            logging.error(f"Unknown kind of variables {port.get_kind()}")
            raise ValueError

        # Set the dimension in the `fem`
        fem.set_dim(dim)

        # Set the mesh of the `domain`, according to the `mesh_ID` of the `port`
        fem.set_mesh(self.domain.get_mesh()[port.get_mesh_id()])

        # Apply the fem to the `port`
        port.set_fem(fem)

        # If `region` is not None, the variable is restricted to the `region` of mesh_id. Useful for boundary ports or interconnected dphs on the same mesh.
        # Add the FEM to the `flow variable` of the `port` in getfem
        if port.get_region() is not None:
            self.gf_model.add_filtered_fem_variable(
                port.get_flow(),
                port.get_fem(),
                port.get_region(),
            )
        else:
            self.gf_model.add_fem_variable(port.get_flow(), port.get_fem())

        # If the port is algebraic AND not substituted, we also add the FEM to the `effort variable` of the `port` in getfem
        if port.get_algebraic() and not port.get_substituted():
            if port.get_region() is not None:
                self.gf_model.add_filtered_fem_variable(
                    port.get_effort(),
                    port.get_fem(),
                    port.get_region(),
                )
            else:
                self.gf_model.add_fem_variable(port.get_effort(), port.get_fem())

        # If the port is dynamic and the co-state is not substituted, the co-state must be added as an unknown variable
        if port.get_effort() in self.costates.keys():
            try:
                assert not port.get_algebraic()
            except AssertionError as err:
                logging.error(
                    f"Port {name_port} should not be algebraic because {port.get_effort()} is a co-state"
                )
                raise err
            if not port.get_substituted():
                if port.get_region() is not None:
                    self.gf_model.add_filtered_fem_variable(
                        port.get_effort(),
                        port.get_fem(),
                        port.get_region(),
                    )
                else:
                    self.gf_model.add_fem_variable(
                        port.get_effort(),
                        port.get_fem(),
                    )

        # If the port is dynamic, the state needs initialization for time-resolution
        if not port.get_algebraic():
            self.initial_value_set[port.get_flow()] = False

    def add_lumped_port(self, name_port: str):
        """Register a global scalar port with no associated finite-element mesh.

        Add the port first, using a State/CoState pair or a Port object with
        kind="lumped". For algebraic ports with substituted=True only the
        flow variable is registered. Dynamic substituted co-states are not
        supported. Global Test_* functions can be used in volumetric GWFL
        forms; divide by the domain measure for an ordinary ODE residual.
        """
        if self.domain is None or not self.domain.get_isSet():
            raise ValueError("Set a domain before adding a lumped port")
        if name_port not in self.ports:
            raise KeyError(f"Unknown port {name_port}")
        port = self.ports[name_port]
        if port.get_kind() != "lumped":
            raise ValueError("Use kind='lumped' for global scalar ports")
        if port.get_region() is not None:
            raise ValueError("A lumped port cannot have a mesh region")
        if port.get_isSet():
            raise ValueError(f"Port {name_port} has already been registered")
        if port.get_substituted() and not port.get_algebraic():
            raise NotImplementedError(
                "Substituted dynamic lumped co-states are unsupported"
            )
        self.gf_model.add_variable(port.get_flow(), 1)
        self.lumped_variables.add(port.get_flow())
        if not port.get_substituted():
            self.gf_model.add_variable(port.get_effort(), 1)
            self.lumped_variables.add(port.get_effort())
        port.set_isSet()
        if not port.get_algebraic():
            self.initial_value_set[port.get_flow()] = False
        if rank == 0:
            logging.info(f"Global scalar port {name_port} registered")

    def add_parameter(self, parameter: Parameter):
        """This function adds a time-independent (possibly space-varying parameter: `x`, `y` and `z` are the space variables) associated to a port of the dphs

        Args:
            parameter (Parameter): the parameter
        """

        name_port = parameter.get_name_port()
        self.ports[name_port].add_parameter(parameter)

        if rank == 0:
            logging.info(
                f"{parameter.get_name()} has been added to port: {parameter.get_name_port()}"
            )

        # Initialization of parameters
        try:
            assert self.ports[name_port].get_isSet()
            self.init_parameter(parameter.get_name(), name_port)
        except AssertionError:
            logging.warning(
                f"{parameter.get_name()} has not been initialized because the port {parameter.get_name_port()} does not have a FEM yet"
            )

    def init_parameter(self, name, name_port):
        """This function initializes the parameter name in the FEM of the port name_port of the dphs and adds it to the getfem `Model`

        Args:
            name (str): the name of the parameter as defined with add_parameter()
            name_port (str): the name of the `port` where the parameter belongs
        """

        # Check if the `port` exists
        try:
            assert name_port in self.ports.keys()
        except AssertionError as err:
            logging.error(f"The port {name_port} does not exist")
            raise err

        # Check if the FEM of the `port` has been set
        try:
            assert self.ports[name_port].get_isSet()
        except AssertionError as err:
            logging.error(
                "The FEM of the port must be set before initializing the parameter"
            )
            raise err

        # Evaluate the expression in the FEM of the `port`
        evaluation = self.ports[name_port].init_parameter(
            name, self.ports[name_port].get_parameter(name).get_expression()
        )

        if rank == 0:
            logging.info(
                f"{name} has been set to {self.ports[name_port].get_parameter(name).get_expression()} in port: {name_port}"
            )

        # Check the size of the parameter according to its kind
        kind = self.ports[name_port].get_parameter(name).get_kind()
        if kind == "tensor-field":
            sizes = evaluation.shape[0]
        elif kind == "scalar-field":
            sizes = 1
        else:
            if rank == 0:
                logging.error(f"Unknown kind of parameter {kind}")
            raise ValueError
        
        # Add the parameter in getfem
        self.gf_model.add_initialized_fem_data(
            name, self.ports[name_port].get_fem(), evaluation, sizes=sizes
        )

        if rank == 0:
            logging.info(
                f"{name} has been initialized with the FEM of port: {name_port}"
            )

    def set_initial_value(
        self, name_variable, expression, *,
        method="l2", projection_mim=None,
    ):
        """Set a field's initial GetFEM coefficients from a Python expression.

        The default `method="l2"` computes a variational L2 projection
        into the actual scalar/vector/derivative FE basis via GetFEM
        assembly and sparse PETSc MUMPS. Historical interpolation remains
        available explicitly via `method="interpolation"`.

        These methods differ for nonpolynomial functions; the default
        is a deliberate v1.3.0 development-branch policy change.
        MPI2, DAE and full physical-model gates remain open.

        Args:
            name_variable (str): initial dynamic variable
            expression (str): SCRIMP Python/NumPy initial expression
            method (str): "l2" (default) or "interpolation" (legacy)
            projection_mim (gf.MeshIm, optional): override the domain's
                integration rule, e.g. HCT composite integration.
        """

        if name_variable in self.initial_value_set.keys():
            if self.initial_value_set[name_variable]:
                logging.warning(f"defining again the initial value of {name_variable}")
        else:
            logging.error(f"Variable {name_variable} does not need initialization")
            raise ValueError

        if method not in ("interpolation", "l2"):
            raise ValueError(
                "Initial FE method must be interpolation or l2"
            )
        if name_variable in self.lumped_variables:
            # Lumped states have no spatial FE basis: preserve their
            # numeric-scalar initialization under either method.
            if name_variable not in self.initial_value_set:
                raise ValueError(f"{name_variable} is not a dynamic state")
            try:
                value = float(expression)
            except (ValueError, TypeError) as exc:
                raise ValueError(
                    "Initial values for lumped states must be numeric scalars"
                ) from exc
            if not np.isfinite(value):
                raise ValueError("Initial value must be finite")
            self.gf_model.set_variable(name_variable, np.array([value]))
            self.initial_value_set[name_variable] = True
            return

        if method == "l2":
            matches = [
                port for port in self.ports.values()
                if name_variable in (port.get_flow(), port.get_effort())
            ]
            if len(matches) != 1:
                raise ValueError(
                    f"Expected exactly one FE port for {name_variable}, "
                    f"found {len(matches)}"
                )
            port = matches[0]
            integration = (
                projection_mim if projection_mim is not None
                else self.domain._int_method[port.get_mesh_id()]
            )
            getfem_rhs_scale = (
                _getfem_mpi_rhs_scale() if max_rank > 1 else 1.0
            )
            values = project_fe_initial_value(
                name_variable, expression, port.get_fem(),
                port.get_region(), integration,
                rhs_scale=getfem_rhs_scale,
            )
            self.set_from_vector(name_variable, values)
            self.initial_value_set[name_variable] = True
            if rank == 0:
                logging.info(
                    f"{name_variable} initialized by GetFEM/PETSc L2 "
                    f"projection from: {expression}"
                )
            return

        evaluation = self.gf_model.mesh_fem_of_variable(name_variable).eval(
            expression, globals(), locals()
        )

        initial_value = None
        for _, port in self.ports.items():
            if name_variable == port.get_flow() or name_variable == port.get_effort():
                if port.get_region() is None:
                    initial_value = evaluation
                else:
                    # A GetFEM MeshFem stores one scalar coefficient for
                    # every intrinsic FE basis function. In particular an
                    # H(div) Raviart–Thomas or BDM vector field does NOT
                    # have qdim independent nodal coefficients: RT0 DOFs
                    # are face-normal fluxes, and their count need not be
                    # divisible by the geometric qdim. Keep the exact
                    # filtered GetFEM coefficient ordering as a 1D array.
                    # The previous reshape (qdim, int(size/qdim)) could
                    # discard a DOF or raise on a genuine RT/BDM basis.
                    fem = port.get_fem()
                    evaluation_long = np.asarray(evaluation).reshape(-1)
                    nb_dofs_total = fem.nbdof()
                    if evaluation_long.size != nb_dofs_total:
                        raise ValueError(
                            f"Initial FEM evaluation of {name_variable} has "
                            f"{evaluation_long.size} coefficients, expected "
                            f"{nb_dofs_total}"
                        )
                    dofs_on_region = np.asarray(
                        fem.basic_dof_on_region(port.get_region())
                    ).reshape(-1)
                    if (not np.issubdtype(dofs_on_region.dtype, np.integer)
                            or np.any(dofs_on_region < 0)
                            or np.any(dofs_on_region >= nb_dofs_total)):
                        raise ValueError(
                            f"Invalid GetFEM region DOF indices for "
                            f"{name_variable}"
                        )
                    initial_value = evaluation_long[dofs_on_region]
                break

        if initial_value is None:
            raise ValueError(f"{name_variable} cannot be found in ports")

        self.set_from_vector(name_variable, initial_value)
        self.initial_value_set[name_variable] = True

        if rank == 0:
            logging.info(f"{name_variable} has been initialized with: {expression}")

    def set_from_vector(self, name_variable, x):
        """This function sets the value of the variable `name_variable` in the getfem `Model` from a numpy vector

        Args:
            name_variable (str): the name of the variable to set
            x (numpy array): the vector of values
        """

        self.gf_model.set_variable(name_variable, x)

        if rank == 0:
            logging.info(f"{name_variable} has been set")
            logging.debug(f"{name_variable} new value is: {x}")

    def add_brick(self, brick: Brick):
        """This function adds a `brick` in the getfem `Model` thanks to a form in GWFL getfem language

        The form may be non-linear

        Args:
            brick (Brick): the brick
        """

        name_brick = brick.get_name()
        mesh_id = brick.get_mesh_id()
        position = brick.get_position()
        form = brick.get_form()

        if brick.get_explicit() and (
            brick.get_dt() or brick.get_linear()
            or brick.get_position() == "source"
        ):
            raise ValueError(
                "explicit=True requires a nonlinear, non-dt, non-source "
                "weak brick; otherwise a lagged RHS cannot be "
                "assembled consistently"
            )
        self.bricks[name_brick] = brick

        # Preserve the original weak-form sign: flow bricks and
        # nonlinear explicit bricks are negated when registered.
        # Explicit terms are subsequently evaluated at the *accepted*
        # TS state and moved to a frozen RHS, never to the Newton tangent.
        if position == "flow" or brick.get_explicit():
            form = "-(" + form + ")"

        for region in brick.get_regions():
            if brick.get_linear() and not brick.get_position() == "source":
                id_brick = self.gf_model.add_linear_term(
                    self.domain._int_method[mesh_id], form, region
                )
                s = "Linear form "

            elif brick.get_position() == "source":
                id_brick = self.gf_model.add_source_term_brick(
                    self.domain._int_method[mesh_id], name_brick[:-7], form, region
                )
                s = "Source form "

            else:
                id_brick = self.gf_model.add_nonlinear_term(
                    self.domain._int_method[mesh_id], form, region
                )
                s = "Non-linear form "

            if rank == 0:
                logging.info(
                    f"{s} '{self.bricks[name_brick].get_form()}' has been added as: {position} relation on region: {region} of mesh: {mesh_id}"
                )

            self.bricks[name_brick].add_id_brick_to_list(id_brick)

            if rank == 0:
                logging.debug(
                    f"brick IDs: {id_brick} has been added to brick: {name_brick}"
                )

    def add_control_port(self, control_port: Control_Port):
        """This function adds a control `port` to the dphs

        Args:
            control_port (Control_Port): the Control Port.
        """

        self.controls[control_port.get_name()] = control_port
        self.add_port(control_port)

    def set_control(self, name, expression):
        """Set a prescribed open-loop input or an implicit feedback law.

        State-independent GWFL expressions use an RHS-only GetFEM source.
        A bare signed model variable (e.g. Y_w or -Y_T) is
        a provably linear feedback: GetFEM assembles its constant
        bilinear coupling directly into the stiffness matrix.
        Other state-dependent expressions, or unknown macros,
        use a nonlinear variational brick so Newton computes
        their full derivative, not merely the residual.

        Args:
            name (str): the name of the control port
            expression (str): GWFL scalar/vector/tensor expression
        """

        kind = self.controls[name].get_kind()
        if kind == "scalar-field":
            times = "*"
        elif kind == "vector-field":
            times = "."
        elif kind == "tensor-field":
            times = ":"
        else:
            logging.error(f"Unknown kind of control {kind}")
            raise ValueError

        # form of the mass matrix for the control variable
        u = self.controls[name].get_name_control()
        mass_form = "-" + u + times + "Test_" + u
        
        self.add_brick(
            Brick(
                "M_" + u,
                mass_form,
                [self.controls[name].get_region()],
                linear=True,
                dt=False,
                position="constitutive",
                mesh_id=self.controls[name].get_mesh_id(),
            )
        )

        # A prescribed control (only t, geometry, known mathematical
        # functions) belongs to the efficient RHS-only source brick.
        #
        # Feedback may depend on another port variable, e.g.
        # set_control("Interface Heat", "Y_w"). A source brick
        # ignores its derivative and would be mathematically WRONG.
        # The exact bare-variable expressions (Y_w, +/-Y_T) are
        # *constant linear* feedback: add_linear_term includes their
        # correct coupling in K and retains the physical Jacobian.
        # All other state-dependent laws stay on add_nonlinear_term
        # with the exact full Newton derivative. Unknown macros are
        # never guessed linear.
        model_variables = {
            variable
            for port in self.ports.values()
            for variable in (port.get_flow(), port.get_effort())
            if variable
        }
        prescribed = _control_expression_is_prescribed(
            expression, model_variables
        )
        direct_linear_feedback = (
            not prescribed
            and _control_expression_is_direct_linear_feedback(
                expression, model_variables
            )
        )
        if prescribed:
            source = Brick(
                u + "_source",
                expression,
                [self.controls[name].get_region()],
                linear=False,
                dt=False,
                position="source",
                mesh_id=self.controls[name].get_mesh_id(),
            )
            # A source generated here is certified state-independent.
            # This tag cannot be inferred from position="source" alone.
            source._scrimp_prescribed_control_source = True
            self.add_brick(source)
        else:
            self.add_brick(
                Brick(
                    u + "_feedback",
                    f"({expression}){times}Test_{u}",
                    [self.controls[name].get_region()],
                    linear=direct_linear_feedback,
                    dt=False,
                    position="constitutive",
                    mesh_id=self.controls[name].get_mesh_id(),
                )
            )

        self.controls[name].get_isSet()
        if rank == 0:
            logging.info(
                f"Control port: {self.controls[name].get_name()} has been set to: {u} = {expression} on region: {self.controls[name].get_region()} of mesh: {self.controls[name].get_mesh_id()}"
            )

    def assemble_mass(self):
        """This function performs the assembly of the bricks dt=True and linear=True and set the PETSc.Mat attribute `mass`"""

        if rank == 0:
            logging.info("Starting linear mass matrix assembly...")

        start = time.perf_counter()

        for _, brick in self.bricks.items():
            # Enable the bricks that are dynamical and linear
            if brick.get_dt() and brick.get_linear():
                brick.enable_id_bricks(self.gf_model)

        _mpi_trace("mass: before GetFEM build_matrix")
        with self.profiler.measure(f"{self._profile_phase}.getfem.mass"):
            self.gf_model.assembly(option="build_matrix")
        _mpi_trace("mass: before GetFEM -> PETSc conversion")
        with self.profiler.measure(f"{self._profile_phase}.bridge.mass"):
            convert_gmm_to_petsc(
                self.gf_model.tangent_matrix(),
                self.mass,
                profiler=self.profiler,
                profile_label=f"{self._profile_phase}.bridge.mass",
                distributed_coo_eligible=True,
            )
        _mpi_trace("mass: PETSc matrix assembled")

        self.disable_all_bricks()

        if rank == 0:
            logging.info(
                f"Linear mass matrix assembly done in {time.perf_counter() - start:1.4g}s"
            )

    def assemble_stiffness(self):
        """This function performs the assembly of the bricks dt=False and linear=True and set the PETSc.Mat attribute `stiffness`"""

        if rank == 0:
            logging.info("Starting linear stiffness matrix assembly...")

        start = time.perf_counter()

        for _, brick in self.bricks.items():
            # Enable the bricks that are non-dynamical and linear
            if not brick.get_dt() and brick.get_linear():
                brick.enable_id_bricks(self.gf_model)

        _mpi_trace("stiffness: before GetFEM build_matrix")
        with self.profiler.measure(f"{self._profile_phase}.getfem.stiffness"):
            self.gf_model.assembly(option="build_matrix")
        _mpi_trace("stiffness: before GetFEM -> PETSc conversion")
        with self.profiler.measure(f"{self._profile_phase}.bridge.stiffness"):
            convert_gmm_to_petsc(
                self.gf_model.tangent_matrix(),
                self.stiffness,
                profiler=self.profiler,
                profile_label=f"{self._profile_phase}.bridge.stiffness",
                distributed_coo_eligible=True,
            )
        _mpi_trace("stiffness: PETSc matrix assembled")

        self.disable_all_bricks()

        if rank == 0:
            logging.info(
                f"Linear stiffness matrix assembly done in {time.perf_counter() - start:1.4g}s"
            )

    def assemble_rhs(self, source_only=False, explicit_only=False):
        """Assemble the GetFEM RHS into the owned PETSc entries.

        GetFEM's model.rhs() mixes two DIFFERENT signs:
        * generic nonlinear weak terms N(z): rhs = -N(z);
        * add_source_term_brick prescribed forcing S(t): rhs = +S(t).

        For raw GetFEM rhs_total=-P*N(z)+S(t), a generic RHS scale
        alpha=1/P alone incorrectly scales the source S/P. The tested
        GetFEM MPI-2 build treats generic terms and source bricks
        differently. Opt-in SCRIMP_GETFEM_SOURCE_RHS_MODE=raw inserts
        unscaled source-only RHS and reconstructs the direct residual as
        -alpha*rhs_total + (1+alpha)*rhs_source = N(z)+S(t).
        match_generic (default) preserves historical uniform scaling.
        This flag only affects source bricks' contribution to the
        direct nonlinear residual, not the physical equations or LU.
        See sandbox/check_mcss_source_sign.py.
        """

        if source_only and explicit_only:
            raise ValueError("RHS assembly cannot be both source-only and explicit-only")
        for _, brick in self.bricks.items():
            if source_only:
                enabled = brick.get_position() == "source"
            elif explicit_only:
                # Freeze the complete nonlinear weak-form contribution
                # at the previous accepted TS state. The registered
                # GetFEM form retains its original sign.
                enabled = brick.get_explicit()
            else:
                # A fully explicit brick must NOT be reevaluated at
                # z_{n+1} in IFunction. Its cached RHS is added later.
                enabled = (
                    brick.get_position() == "source"
                    or (
                        not brick.get_dt()
                        and not brick.get_linear()
                        and not brick.get_explicit()
                    )
                )
            if enabled:
                brick.enable_id_bricks(self.gf_model)

        _mpi_trace("rhs: before GetFEM build_rhs")
        with self.profiler.measure(f"{self._profile_phase}.getfem.rhs"):
            self.gf_model.assembly(option="build_rhs")
        _mpi_trace("rhs: before PETSc owned-row insertion")
        # GetFEM MPI-2 already sums and replicates the entire RHS on
        # every rank. Insert only PETSc-owned rows; never ADD these
        # replicated values, which would multiply them by the rank count.
        gf_rhs = np.asarray(self.gf_model.rhs())
        size = self.gf_model.nbdof()
        if gf_rhs.size != size or self.rhs.getSize() != size:
            raise ValueError("GetFEM/PETSc RHS size mismatch")
        # For the local GetFEM 5.4.2 MPI-2 installation, generic
        # terms' *global* RHS is P-fold too large at P=1,2,4;
        # the native matrices are not. This indicates duplicate
        # contributions before/in GetFEM's MPI reduction, not duplicate
        # PETSc insertion (only owned rows are copied here).
        #
        # WARNING: full RHS replication is normal and is NOT itself
        # a reason to divide by P. 'auto' runs a small exact-integral
        # calibration once per MPI process and scales *only if* a
        # duplicate factor is observed. 'raw' is for already-correct
        # GetFEM RHS or model-specific exceptions. Explicit RHS/
        # constraint bricks have special master-only GetFEM semantics:
        # a generic probe does not certify every such mixed model.
        rhs_mode = os.environ.get(
            "SCRIMP_GETFEM_RHS_MODE", "auto"
        ).strip().lower()
        if rhs_mode not in ("auto", "raw", "divide_by_mpi_size"):
            raise ValueError(
                "SCRIMP_GETFEM_RHS_MODE must be "
                "auto, raw or divide_by_mpi_size"
            )
        if rhs_mode == "auto":
            scale = _getfem_mpi_rhs_scale()
        elif rhs_mode == "divide_by_mpi_size":
            scale = 1.0 / max_rank
        else:
            scale = 1.0
        source_mode = _getfem_source_rhs_mode()
        if not source_only and not explicit_only:
            # Save alpha used for the full mixed (-N+S) RHS. The direct
            # residual needs it to undo source attenuation when source
            # bricks are already globally assembled under GetFEM MPI-2.
            self._getfem_full_rhs_scale = scale
        elif source_only and source_mode == "raw":
            # On the validated GetFEM MPI-2 build, SOURCE bricks are
            # globally summed ONCE, unlike generic nonlinear terms.
            # The separately frozen explicit nonlinear RHS is still
            # a generic term and MUST retain the calibrated scale.
            scale = 1.0
        if os.environ.get("SCRIMP_MPI_TRACE") == "1":
            _mpi_trace(
                f"rhs: mode={rhs_mode}, source_mode={source_mode}, "
                f"source_only={source_only}, scale={scale:.9g}"
            )
        with self.profiler.measure(f"{self._profile_phase}.bridge.rhs"):
            lo, hi = self.rhs.getOwnershipRange()
            if hi > lo:
                self.rhs.setValues(
                    range(lo, hi),
                    gf_rhs[lo:hi] * scale,
                    addv=PETSc.InsertMode.INSERT_VALUES,
                )
            self.rhs.assemble()
        _mpi_trace("rhs: PETSc vector assembled")

        self.disable_all_bricks()

    def _freeze_explicit_rhs(self, time_at_accepted_state, accepted_z):
        """Freeze all fully explicit nonlinear bricks at one known TS state.

        This is a first-order lagged RHS: E(t_n, z_n), independent of
        Newton trial state z_{n+1} and trial derivative. Call this once
        at the *beginning* of each accepted-state step, including the
        historical pseudo/BDF initialization TS. A rejected attempt
        must reuse the same accepted-state forcing.
        """
        if not any(brick.get_explicit() for brick in self.bricks.values()):
            return
        if self._explicit_rhs is None:
            self._explicit_rhs = self.rhs.duplicate()
        self._explicit_rhs_ready = False
        self._sync_lumped_trial_state(time_at_accepted_state, accepted_z)
        self.assemble_rhs(explicit_only=True)
        self.rhs.copy(self._explicit_rhs)
        self._explicit_rhs_ready = True
        self._explicit_rhs_time = float(time_at_accepted_state)

    def _prestep_explicit_rhs(self, ts):
        """PETSc TS PreStep: use ONLY the last accepted state, not SNES."""
        if not any(brick.get_explicit() for brick in self.bricks.values()):
            return
        self._freeze_explicit_rhs(ts.getTime(), ts.getSolution())

    def assemble_nl_mass(self):
        """This function performs the assembly of the bricks dt=True and linear=False and set the PETSc.Mat attribute `nl_mass`"""

        if self._matrix_memory_mode == "compact" and self.linear_mass:
            # The tangent aliases the *assembled* linear mass; do not
            # accidentally change that read-only matrix via zeroEntries.
            return

        for _, brick in self.bricks.items():
            # Enable the bricks that are dynamical and non-linear
            if brick.get_dt() and not brick.get_linear():
                brick.enable_id_bricks(self.gf_model)

        _mpi_trace("nl_mass: before GetFEM build_matrix")
        with self.profiler.measure(f"{self._profile_phase}.getfem.nl_mass"):
            self.gf_model.assembly(option="build_matrix")
        _mpi_trace("nl_mass: before GetFEM -> PETSc conversion")
        with self.profiler.measure(f"{self._profile_phase}.bridge.nl_mass"):
            convert_gmm_to_petsc(
                self.gf_model.tangent_matrix(),
                self.nl_mass,
                profiler=self.profiler,
                profile_label=f"{self._profile_phase}.bridge.nl_mass",
                distributed_coo_eligible=True,
            )
        _mpi_trace("nl_mass: PETSc matrix assembled")
        _mpi_trace("nl_mass: before PETSc copy/AXPY")
        # Keep the allocated PETSc Mat alive throughout TS/Newton iterations.
        # Reassigning self.tangent_mass = self.mass.copy() at every residual
        # evaluation creates many temporary MPIAIJ matrices. In a 2-rank
        # shallow-water run these were freed by the TS monitor's forced
        # PetscGarbageCleanup, triggering an invalid-free in MatDestroy_MPIAIJ.
        # Reassemble in place; the nonlinear sparsity may change with z.
        structure = PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN
        self.tangent_mass.setOption(
            PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False
        )
        try:
            self.tangent_mass.zeroEntries()
            self.tangent_mass.axpy(1, self.mass, structure=structure)
            _mpi_trace("nl_mass: tangent reused; before nonlinear AXPY")
            self.tangent_mass.axpy(1, self.nl_mass, structure=structure)
            self.tangent_mass.assemble()
        finally:
            self.tangent_mass.setOption(
                PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, True
            )
        _mpi_trace("nl_mass: tangent reassembled in place")

        self.disable_all_bricks()
        _mpi_trace("nl_mass: before barrier")
        comm.barrier()
        _mpi_trace("nl_mass: after barrier")

    def assemble_nl_stiffness(self):
        """This function performs the assembly of the bricks dt=False and linear=False and set the PETSc.Mat attribute `nl_stiffness`"""

        if self._matrix_memory_mode == "compact" and self.linear_stiffness:
            # Do not modify the aliased linear stiffness by accidentally
            # assembling an empty nonlinear tangent.
            return

        for _, brick in self.bricks.items():
            # Enable the bricks that are non-dynamical and non-linear
            if (
                not brick.get_dt()
                and not brick.get_linear()
                and not brick.get_explicit()
            ):
                brick.enable_id_bricks(self.gf_model)

        _mpi_trace("nl_stiffness: before GetFEM build_matrix")
        with self.profiler.measure(f"{self._profile_phase}.getfem.nl_stiffness"):
            self.gf_model.assembly(option="build_matrix")
        _mpi_trace("nl_stiffness: before GetFEM -> PETSc conversion")
        with self.profiler.measure(f"{self._profile_phase}.bridge.nl_stiffness"):
            convert_gmm_to_petsc(
                self.gf_model.tangent_matrix(),
                self.nl_stiffness,
                profiler=self.profiler,
                profile_label=f"{self._profile_phase}.bridge.nl_stiffness",
                distributed_coo_eligible=True,
            )
        _mpi_trace("nl_stiffness: PETSc matrix assembled")
        _mpi_trace("nl_stiffness: before PETSc copy/AXPY")
        # As for the nonlinear mass, retain the existing PETSc Mat handle.
        # Newton iteration must not enqueue an obsolete MPIAIJ matrix for
        # garbage collection at every assembly.
        structure = PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN
        self.tangent_stiffness.setOption(
            PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False
        )
        try:
            self.tangent_stiffness.zeroEntries()
            self.tangent_stiffness.axpy(1, self.stiffness, structure=structure)
            _mpi_trace("nl_stiffness: tangent reused; before nonlinear AXPY")
            self.tangent_stiffness.axpy(
                1, self.nl_stiffness, structure=structure
            )
            self.tangent_stiffness.assemble()
        finally:
            self.tangent_stiffness.setOption(
                PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, True
            )
        _mpi_trace("nl_stiffness: tangent reassembled in place")

        self.disable_all_bricks()
        _mpi_trace("nl_stiffness: before barrier")
        comm.barrier()
        _mpi_trace("nl_stiffness: after barrier")

    def disable_all_bricks(self):
        """This function disables all bricks in the `Model`"""
        for _, brick in self.bricks.items():
            brick.disable_id_bricks(self.gf_model)

    def enable_all_bricks(self):
        """This function enables all bricks in the `Model`"""
        for _, brick in self.bricks.items():
            brick.enable_id_bricks(self.gf_model)

    def _sync_lumped_trial_state(self, t, z):
        """Synchronize GetFEM's time and nonlinear Newton trial state.

        Despite the historical name, nonlinear *distributed* PDEs also
        require the current PETSc iterate (not the previous TS monitor state).
        All processes enter the collective scatter in the nonlinear case.
        """
        self.gf_model.set_time(float(t))
        if self.linear_mass and self.linear_stiffness and not self.lumped_variables:
            return
        self.Scatter[0].scatter(z, self.GlobalVec)
        # PETSc is compiled with complex scalars, while GetFEM's
        # physical model variables are real. Select the real component
        # explicitly, as in the TS monitor, instead of letting NumPy
        # silently cast complex128 -> float64 (ComplexWarning).
        # Preserve the historical copy semantics of the cast.
        self.gf_model.to_variables(
            np.array(
                self.GlobalVec.getArray(readonly=True).real,
                dtype=np.float64,
                copy=True,
            )
        )

    def _trace_newton_blocks(self, phase, t, z, zd):
        """Optional Newton diagnostics by physical/GetFEM variable.

        GetFEM exposes a complete global trial vector after our existing
        VecScatter. This is read-only and rank-0-only: no new collectives,
        new matrix assembly or changes to the Newton iteration.
        """
        if os.environ.get("SCRIMP_NEWTON_BLOCK_TRACE") != "1" or rank != 0:
            return
        state = np.asarray(
            self.GlobalVec.getArray(readonly=True)
        ).real
        # GetFEM Model.variable_list() is a print-only command in
        # GetFEM 5.4.2; it returns None, so it cannot be iterated.
        # SCRIMP's ports are our authoritative registry of all solver
        # unknowns (flow and, when present, effort). Query GetFEM only
        # for the global index interval of each registered variable.
        # dict.fromkeys preserves registration order and removes any
        # duplicate names (e.g. a substituted flow/effort pair).
        names = dict.fromkeys(
            name
            for port in self.ports.values()
            for name in (port.get_flow(), port.get_effort())
            if name
        )
        rows = []
        for name in names:
            try:
                start, width = self.gf_model.interval_of_variable(name)
            except (KeyError, ValueError, RuntimeError):
                # A substituted effort or non-unknown is not present
                # in the GetFEM solver vector. Do not try to print it.
                continue
            block = state[int(start):int(start + width)]
            if not block.size:
                continue
            norm = float(np.linalg.norm(block))
            max_abs = float(np.max(np.abs(block)))
            annotation = ""
            if name == "h":
                annotation = (
                    f", min={np.min(block):.6g}, "
                    f"max={np.max(block):.6g}"
                )
            rows.append(
                f"{name}:norm={norm:.6g},maxabs={max_abs:.6g}"
                f"{annotation}"
            )
        # These norms are PETSc collectives; run them on *all* ranks
        # elsewhere if needed. Here the globally replicated GetFEM
        # state supplies all physically meaningful per-variable norms.
        print(
            f"[SCRIMP Newton blocks] {phase} t={t:.9g} "
            f"global_norm={np.linalg.norm(state):.6g} "
            + " | ".join(rows),
            file=sys.stderr, flush=True,
        )

    def IFunction(self, TS, t, z, zd, F):
        """
        IFunction for the time-resolution of the dphs with PETSc TS fully implicit integrator

        Args:
            TS (PETSc TS): the PETSc TS object handling time-resolution
            t (float): time parameter
            z (PETSc Vec): the state
            zd (PETSc Vec): the time-derivative of the state
            F (PETSc Vec): the rhs vector
        """

        _mpi_trace(f"IFunction: enter t={t:.8g}")
        self._sync_lumped_trial_state(t, z)
        _mpi_trace("IFunction: Newton state synchronized")
        self._trace_newton_blocks("IFunction", t, z, zd)
        if not self.linear_mass:
            self.assemble_nl_mass()
        nonlinear_form = os.environ.get(
            "SCRIMP_NONLINEAR_RESIDUAL_FORM", "direct_getfem"
        ).strip().lower()
        # The direct weak-form residual does NOT need the nonlinear
        # tangent matrix. Avoid expensive GetFEM/PETSc sparse matrix
        # conversions at every function evaluation (including line
        # search trials); IJacobian assembles the tangent when needed.
        if not self.linear_stiffness and nonlinear_form == "legacy":
            self.assemble_nl_stiffness()

        if any(brick.get_explicit() for brick in self.bricks.values()):
            if not self._explicit_rhs_ready or self._explicit_rhs is None:
                raise RuntimeError(
                    "explicit=True nonlinear bricks require a frozen RHS "
                    "from the previous accepted state; use PETSc TS "
                    "PreStep or _freeze_explicit_rhs() in a standalone "
                    "Jacobian diagnostic"
                )
        self.assemble_rhs()
        _mpi_trace("IFunction: RHS assembled; before PETSc algebra")

        if nonlinear_form not in ("legacy", "direct_getfem"):
            raise ValueError(
                "SCRIMP_NONLINEAR_RESIDUAL_FORM must be "
                "legacy or direct_getfem"
            )

        if nonlinear_form == "direct_getfem":
            # Correct general residual: GetFEM build_rhs is the NEGATIVE
            # of the nonlinear weak-form residual. K_nl(z) is its
            # derivative, NOT a nonlinear matrix to multiply by z.
            # Validated on full viscous Saint-Venant 2D, 1/2 MPI ranks,
            # for independent z and zd perturbations.
            mass_matrix = (
                self.mass if self.linear_mass else self.tangent_mass
            )
            mass_matrix.mult(zd, self.buffer)
            self.stiffness.multAdd(z, self.buffer, F)
            F.axpy(-1, self.rhs)
            if any(
                brick.get_position() == "source"
                for brick in self.bricks.values()
            ):
                # The full scaled RHS is alpha*(-P*N+S); with source
                # raw=+S (opt-in), -full + (1+alpha)*source = N+S.
                # Without opt-in, source is also scaled by alpha:
                # -full + 2*source preserves the previous behavior.
                source_mode = _getfem_source_rhs_mode()
                correction = (
                    1.0 + self._getfem_full_rhs_scale
                    if source_mode == "raw" else 2.0
                )
                self.assemble_rhs(source_only=True)
                F.axpy(correction, self.rhs)
        elif self.lumped_variables and self.linear_mass:
            # GetFEM's build_rhs returns MINUS the nonlinear weak-form
            # residual (plus source contributions).  Its nonlinear tangent
            # matrix is the Jacobian, not an operator to multiply by the
            # current trial state when evaluating the residual.
            #
            # The legacy K_tangent*z + rhs reconstruction only happens to
            # be correct for certain homogeneous expressions (notably
            # quadratic ones).  It reverses the sign of constant forcing
            # and degree-zero rational constitutive laws, as exposed by
            # the fluid/tank PDE-ODE coupling.  Evaluate the nonlinear
            # weak form using -rhs, while keeping K_tangent for IJacobian.
            self.mass.mult(zd, self.buffer)
            self.stiffness.multAdd(z, self.buffer, F)
            F.axpy(-1, self.rhs)
            if any(
                brick.get_position() == "source"
                for brick in self.bricks.values()
            ):
                source_mode = _getfem_source_rhs_mode()
                correction = (
                    1.0 + self._getfem_full_rhs_scale
                    if source_mode == "raw" else 2.0
                )
                self.assemble_rhs(source_only=True)
                F.axpy(correction, self.rhs)
        else:
            # Legacy reproduction for historical scripts; use only
            # to investigate older results, not for correct Newton.
            self.tangent_mass.mult(zd, self.buffer)
            self.tangent_stiffness.multAdd(z, self.buffer, F)
            F.axpy(1, self.rhs)
        if self._explicit_rhs_ready:
            if nonlinear_form != "direct_getfem":
                raise ValueError(
                    "Fully explicit bricks require "
                    "SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem"
                )
            # For a nonlinear weak brick GetFEM RHS=-N. The frozen
            # value has the SAME sign as the original IFunction,
            # but is independent of the current Newton iteration.
            F.axpy(-1, self._explicit_rhs)
        _mpi_trace("IFunction: residual completed")

    def _add_nonlinear_mass_hessian_fd(self, t, z, zd, P):
        """Add (D_z M(z)[.] ) zd to PETSc's implicit Jacobian.

        Nonlinear GetFEM dt bricks represent a weak flow primitive
        Q(z), with M(z) = D Q(z). By Schwarz symmetry,
        (D_z M(z)[v]) zd = (D_z M(z)[zd]) v, which can be assembled
        by a directional central difference of the tangent mass.

        Default correction for nonlinear time bricks; its symmetric
        finite-difference approximation was verified on the full
        2D viscous Saint-Venant Jacobian to ~1e-11. The legacy option
        remains available for compatibility comparisons.
        """
        if self.linear_mass or zd.norm() == 0:
            return

        eta = float(os.environ.get(
            "SCRIMP_NONLINEAR_MASS_FD_EPSILON", "1e-4"
        ))
        if not np.isfinite(eta) or eta <= 0:
            raise ValueError(
                "SCRIMP_NONLINEAR_MASS_FD_EPSILON must be finite > 0"
            )

        # The TS pseudo-initialization may use dt ~ 1e-9, so its
        # derivative zd=(z_new-z_old)/dt can be orders of magnitude
        # larger than the physical state. The raw z +/- 1e-4*zd
        # probe must NOT jump to enormous, nonphysical heights.
        # Keep ||eta*zd|| <= epsilon * max(1,||z||); for ordinary
        # moderate zd this retains the already FD-validated eta.
        zd_norm = float(zd.norm())
        z_norm = float(z.norm())
        if not np.isfinite(zd_norm) or not np.isfinite(z_norm):
            raise FloatingPointError(
                "Non-finite state or time-derivative norm in "
                "nonlinear mass Hessian assembly"
            )
        eta *= min(
            1.0,
            max(1.0, z_norm) / max(1.0, zd_norm),
        )
        _mpi_trace(
            "mass_hessian: z_norm={:.6g}, zd_norm={:.6g}, "
            "eta={:.6g}, perturb_norm={:.6g}".format(
                z_norm, zd_norm, eta, eta * zd_norm,
            )
        )

        if self._matrix_memory_mode == "compact":
            # One reusable PETSc Mat and Vec across IJacobian calls.
            # The - perturbation remains in tangent_mass; no extra
            # copy or correction Mat is necessary.
            add_mass_hessian_fd_reused(
                self, t, z, zd, P, eta,
                PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN,
            )
            return

        # Original allocation-heavy path retained for A/B parity.
        z_plus, z_minus = z.copy(), z.copy()
        mass_plus = mass_minus = correction = None
        try:
            z_plus.axpy(eta, zd)
            z_minus.axpy(-eta, zd)

            self._sync_lumped_trial_state(t, z_plus)
            self.assemble_nl_mass()
            mass_plus = self.tangent_mass.copy()

            self._sync_lumped_trial_state(t, z_minus)
            self.assemble_nl_mass()
            mass_minus = self.tangent_mass.copy()

            structure = PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN
            correction = mass_plus.copy()
            correction.axpy(-1.0, mass_minus, structure=structure)
            correction.scale(0.5 / eta)
            P.axpy(1.0, correction, structure=structure)
            P.assemble()
        finally:
            # Crucial: callbacks after IJacobian must see the original
            # Newton trial state, not the displaced Hessian state.
            self._sync_lumped_trial_state(t, z)
            self.assemble_nl_mass()
            for obj in (correction, mass_plus, mass_minus,
                        z_plus, z_minus):
                if obj is not None:
                    obj.destroy()

    def set_linear_flags(self):
        """Flag state-dependent mass/stiffness terms for reassembly.

        GetFEM source bricks contribute only to the RHS, even when
        declared linear=False for their special assembly route.
        Feedback controls created by set_control() are *not* sources:
        signed bare-variable feedback is a constant linear coupling
        assembled in K; general feedback has a state-dependent
        constitutive tangent and is reassembled per Newton iterate.
        """
        self.linear_mass = True
        self.linear_stiffness = True
        for brick in self.bricks.values():
            if brick.get_position() == "source":
                continue
            if brick.get_dt() and not brick.get_linear():
                self.linear_mass = False
            if not brick.get_dt() and not brick.get_linear():
                self.linear_stiffness = False

    def IJacobian(self, TS, t, z, zd, sig, A, P):
        """
        IJacobian for the time-resolution of the dphs with PETSc TS fully implicit integrator

        Args:
            TS (PETSc TS): the PETSc TS object handling time-resolution
            t (float): time parameter
            z (PETSc Vec): the state
            zd (PETSc Vec): the time-derivative of the state
            sig (float): a shift-parameter, depends on dt
            A (PETSc Mat): the jacobian matrix
            P:(PETSc Mat) the jacobian matrix to use for pre-conditionning
        """

        _mpi_trace(f"IJacobian: enter t={t:.8g}, shift={sig:.8g}")
        self._sync_lumped_trial_state(t, z)
        _mpi_trace("IJacobian: Newton state synchronized")
        self._trace_newton_blocks("IJacobian", t, z, zd)
        if not self.linear_mass:
            self.assemble_nl_mass()
        if not self.linear_stiffness:
            self.assemble_nl_stiffness()

        # In a nonlinear PDE the sparsity pattern can change between
        # Newton evaluations. P is created from the *initial* tangent
        # and may therefore have a smaller pattern. MatCopy with
        # SAME_NONZERO_PATTERN on matrices with different sparsity is
        # explicitly unsafe in PETSc (and can crash). Likewise MatAXPY
        # must not assume that M(z) has the same nonzeros as K(z).
        structure = PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN
        # v1.3.0 default: on_change tracks the ACTUAL assembled PETSc
        # operator, not model linearity. Reuse the symbolic LU when its
        # exact stored graph is unchanged, and the numeric LU only when
        # all stored coefficients are identical across all MPI ranks.
        # This handles nonlinear tangents and changing TS shift correctly.
        # Explicit never/always policies remain supported as overrides.
        # Checks happen after ALL Jacobian contributions are assembled.
        policy = getattr(self, "_mumps_reset_policy", None)
        if policy is None:
            policy = _get_mumps_reset_policy(
                all_bricks_linear=_all_bricks_are_linear(self.bricks)
            )
        ksp = TS.getSNES().getKSP()
        pc = ksp.getPC()
        is_mumps = (
            pc.getType() == "lu"
            and pc.getFactorSolverType() == "mumps"
        )
        if policy == "always" and is_mumps:
            _mpi_trace("IJacobian: always policy before MUMPS PC.reset")
            with self.profiler.measure(
                f"{self._profile_phase}.petsc.pc_reset"
            ):
                pc.reset()
            _mpi_trace("IJacobian: always policy after MUMPS PC.reset")

        # SNESComputeJacobian() overwrites KSP's reuse-preconditioner
        # flag AFTER this callback returns (PETSc 3.26). Merely calling
        # pc.setReusePreconditioner(True) here is therefore ineffective.
        # Instead, assemble into a reusable *independent* workspace.
        # If the full CSR operator is unchanged, do NOT touch PETSc's
        # factorized P (or A): its MatState remains unchanged, allowing
        # PCSetUp() to skip numerical LU without any SNES lagging.
        # In particular, this also works if the TS shift or nonlinear
        # tangent changes: a changed operator is copied back.
        compare_on_change = policy == "on_change" and is_mumps
        if compare_on_change:
            candidate = getattr(self, "_mumps_jacobian_candidate", None)
            if candidate is None:
                # No SHARED_NONZERO_PATTERN: nonlinear operators may gain
                # or lose entries and must remain independently mutable.
                candidate = self.tangent_stiffness.duplicate(copy=True)
                self._mumps_jacobian_candidate = candidate
            target = candidate
        else:
            target = P

        _mpi_trace("IJacobian: before tangent -> Jacobian work matrix")
        self.tangent_stiffness.copy(target, structure=structure)
        _mpi_trace("IJacobian: after tangent -> Jacobian work matrix")
        _mpi_trace("IJacobian: before shift*mass AXPY")
        target.axpy(sig, self.tangent_mass, structure=structure)

        mass_mode = os.environ.get(
            "SCRIMP_NONLINEAR_MASS_JACOBIAN", "finite_difference"
        ).strip().lower()
        if mass_mode not in ("legacy", "finite_difference"):
            raise ValueError(
                "SCRIMP_NONLINEAR_MASS_JACOBIAN must be "
                "legacy or finite_difference"
            )
        if mass_mode == "finite_difference" and not self.linear_mass:
            if os.environ.get(
                "SCRIMP_NONLINEAR_RESIDUAL_FORM", "direct_getfem"
            ).strip().lower() != "direct_getfem":
                raise ValueError(
                    "The mass Hessian diagnostic requires "
                    "SCRIMP_NONLINEAR_RESIDUAL_FORM=direct_getfem"
                )
            _mpi_trace("IJacobian: before nonlinear mass Hessian FD")
            self._add_nonlinear_mass_hessian_fd(t, z, zd, target)
            _mpi_trace("IJacobian: after nonlinear mass Hessian FD")

        _mpi_trace("IJacobian: preconditioner assembled")

        if os.environ.get("SCRIMP_MPI_TRACE") == "1":
            # Local-only diagnostics, with no extra MPI collectives.
            # This catches singular/non-finite operators before PETSc's
            # native KSP/PC factorization path is entered.
            rp, ci, vals = target.getValuesCSR()
            bad = int(np.count_nonzero(~np.isfinite(vals)))
            diagonal = target.getDiagonal()
            d = np.asarray(diagonal.getArray(readonly=True))
            diag_zeros = int(np.count_nonzero(d == 0))
            diag_bad = int(np.count_nonzero(~np.isfinite(d)))
            # Zero diagonals may be expected in mixed saddle-point DAEs.
            # Entirely zero numerical rows, however, indicate an absent
            # equation. Detect the latter separately (trace mode only).
            first_row = target.getOwnershipRange()[0]
            numerical_zero_rows = np.asarray([
                not np.any(vals[rp[i]:rp[i + 1]])
                for i in range(len(rp) - 1)
            ], dtype=bool)
            zero_rows_global = (
                np.flatnonzero(numerical_zero_rows) + first_row
            )
            zero_diagonal_global = np.flatnonzero(d == 0) + first_row
            _mpi_trace(
                f"IJacobian: operator shape={target.getSize()}, "
                f"ownership={target.getOwnershipRange()}, "
                f"nnz_owned={len(vals)}, nonfinite={bad}, "
                f"diagonal_zeros={diag_zeros}, diagonal_nonfinite={diag_bad}, "
                f"zero_numerical_rows={len(zero_rows_global)}; "
                f"zero_row_ids={zero_rows_global[:16].tolist()}, "
                f"zero_diag_ids={zero_diagonal_global[:16].tolist()}"
            )
            diagonal.destroy()
            if bad or diag_bad:
                raise FloatingPointError("Non-finite PETSc Jacobian entries")

        if not compare_on_change and A != P:
            # Preserve the original behavior for non-MUMPS/always/never.
            _mpi_trace("IJacobian: copy preconditioner into distinct operator")
            P.copy(A, structure=structure)

        if policy == "on_change":
            # Initial and main TS own distinct KSP/PC/MUMPS factors.
            # Check solver uniformity once for each phase, then do one
            # collective bitmask per Jacobian (1=graph, 2=coefficients).
            mpi = comm.tompi4py()
            phase = self._profile_phase
            checked = getattr(self, "_mumps_phase_enabled", None)
            if checked is None:
                checked = self._mumps_phase_enabled = {}
            if phase not in checked:
                n_mumps = mpi.allreduce(int(is_mumps), op=MPI.SUM)
                if n_mumps not in (0, max_rank):
                    raise RuntimeError(
                        "MUMPS LU must be configured on every MPI rank "
                        "for collective on_change symbolic/numeric reuse"
                    )
                checked[phase] = n_mumps == max_rank
            if is_mumps != checked[phase]:
                raise RuntimeError(
                    "MUMPS solver type changed within the same TS phase; "
                    "this invalidates the collective reuse cache"
                )
            if is_mumps:
                with self.profiler.measure(
                    f"{self._profile_phase}.petsc.csr_pattern_check"
                ):
                    try:
                        current = _local_stored_csr_snapshot(target)
                        csr_error = None
                    except (PETSc.Error, AttributeError,
                            NotImplementedError, TypeError) as exc:
                        current = None
                        csr_error = exc
                    if mpi.allreduce(
                        int(csr_error is not None), op=MPI.SUM
                    ):
                        raise RuntimeError(
                            "Cannot inspect full stored PETSc CSR on "
                            "every rank; automatic symbolic/numeric "
                            "reuse cannot be validated"
                        ) from csr_error
                    snapshots = getattr(self, "_mumps_csr_snapshots", None)
                    if snapshots is None:
                        snapshots = self._mumps_csr_snapshots = {}
                    previous = snapshots.get(phase)
                    local_change = _local_csr_change_flags(previous, current)
                    # MPI BOR broadcasts each kind of change to ALL ranks.
                    global_change = int(
                        mpi.allreduce(local_change, op=MPI.BOR)
                    )
                    graph_changed = bool(global_change & 1)
                    values_changed = bool(global_change & 2)
                    if previous is None or values_changed:
                        snapshots[phase] = current
                # Lightweight per-phase accounting for reproducible
                # qualification. These are collective decisions derived
                # from global_change, not rank-local guesses. Retain
                # only counters after solve, not any CSR data.
                counters = getattr(self, "_mumps_reuse_stats", None)
                if counters is None:
                    counters = self._mumps_reuse_stats = {}
                counts = counters.setdefault(phase, {
                    "ijacobian_calls": 0,
                    "numeric_updates": 0,
                    "symbolic_updates": 0,
                    "reused_operators": 0,
                })
                counts["ijacobian_calls"] += 1
                if previous is None or values_changed:
                    counts["numeric_updates"] += 1
                    if previous is None or graph_changed:
                        counts["symbolic_updates"] += 1
                else:
                    counts["reused_operators"] += 1
                if os.environ.get("SCRIMP_MUMPS_REUSE_DIAG") == "1":
                    # Diagnostic ONLY: no tolerance-based reuse decision.
                    # PETSc CN can change its shift, including on the
                    # final MATCHSTEP interval. We need to distinguish
                    # a genuine matrix change from a redundant PC setup.
                    # One extra collective is only active for this
                    # short opt-in qualification.
                    if (previous is not None
                            and not (local_change & 1)
                            and current[3].shape == previous[3].shape):
                        delta = np.abs(current[3] - previous[3])
                        local_max_diff = (
                            float(np.max(delta)) if delta.size else 0.0
                        )
                        local_different = int(np.count_nonzero(
                            current[3] != previous[3]
                        ))
                    else:
                        local_max_diff = -1.0
                        local_different = -1
                    global_max_diff = float(
                        mpi.allreduce(local_max_diff, op=MPI.MAX)
                    )
                    global_different = int(
                        mpi.allreduce(max(0, local_different), op=MPI.SUM)
                    )
                    if rank == 0:
                        print(
                            "MUMPS-REUSE-AUDIT "
                            f"phase={phase} "
                            f"step={TS.getStepNumber()} "
                            f"time={t:.17g} shift={sig:.17g} "
                            f"graph_changed={int(graph_changed)} "
                            f"values_changed={int(values_changed)} "
                            f"local_flags={local_change} "
                            f"global_flags={global_change} "
                            f"entries_changed_global={global_different} "
                            f"max_abs_delta_global={global_max_diff:.17g} "
                            f"update_operator={int(previous is None or values_changed)}",
                            flush=True,
                        )
                if previous is None or values_changed:
                    # Only a genuinely changed operator may increment the
                    # MatState of P/A. PETSc tests MatState in PCSetUp()
                    # and retains the LU if it has not changed.
                    if previous is not None and graph_changed:
                        _mpi_trace(
                            "IJacobian: stored CSR graph changed globally; "
                            "resetting MUMPS symbolic factor"
                        )
                        with self.profiler.measure(
                            f"{self._profile_phase}.petsc.pc_reset"
                        ):
                            pc.reset()
                    target.copy(P, structure=structure)
                    if A != P:
                        P.copy(A, structure=structure)
                    if previous is not None and graph_changed:
                        ksp.setOperators(A, P)
                # Otherwise DO NOTHING to A or P, including no assembly,
                # no MatCopy, no MatAXPY, no forced PC flag. SNES may
                # explicitly set reuse_preconditioner=False and PCSetUp
                # still skips setup for an identical MatState.

        # Opt-in native diagnostic: dump selected actual SWE Jacobians
        # before PETSc/MUMPS factorization. Binary MatView is collective,
        # so EVERY MPI rank takes part in the same call. Persist the
        # *assembled* preconditioner/operator used by LU. This allows
        # replay in a different PETSc/MUMPS environment without GetFEM.
        #
        # Power-of-two snapshots bound disk I/O during divergent Newton:
        # 1, 2, 4, 8, 16, 32 (up to six files). Normally disabled.
        snapshot_prefix = os.environ.get(
            "SCRIMP_JACOBIAN_SNAPSHOT_PREFIX", ""
        ).strip()
        if snapshot_prefix:
            count = getattr(self, "_jacobian_snapshot_count", 0) + 1
            self._jacobian_snapshot_count = count
            if count in (1, 2, 4, 8, 16, 32):
                directory = os.path.dirname(snapshot_prefix)
                if rank == 0 and directory:
                    os.makedirs(directory, exist_ok=True)
                comm.barrier()
                name = f"{snapshot_prefix}_{count:03d}.petscbin"
                _mpi_trace(f"IJacobian: dumping PETSc binary snapshot {name}")
                viewer = PETSc.Viewer().createBinary(
                    name, mode="w", comm=comm
                )
                try:
                    P.view(viewer)
                finally:
                    viewer.destroy()
                _mpi_trace(f"IJacobian: snapshot {count} written")
        _mpi_trace("IJacobian: completed")

    def init_matrix(self,A):
        """
        Set all diagonal entries of A without modifying existing entries.

        Parameters
        ----------
        A : PETSc.Mat

        """

        A.setSizes(self.gf_model.nbdof())
        A.setOption(PETSc.Mat.Option.FORCE_DIAGONAL_ENTRIES, True)
        A.setType("aij")
        A.setUp()
        # Vec.setArray() expects an array of local values, not scalar 0.
        # Vec.set(0.0) initializes the existing distributed vector safely.
        x = A.createVecLeft()
        x.set(0.0)
        A.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False)
        try:
            A.setDiagonal(x, addv=PETSc.InsertMode.ADD_VALUES)
            A.assemble()
        finally:
            A.setOption(PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, True)
            x.destroy()

    def allocate_memory(self):
        """Pre-allocate memory for matrices and vectors"""

        self.mass.bindToCPU(True)
        self.init_matrix(self.mass)
        self.init_matrix(self.nl_mass)
        if self._matrix_memory_mode == "compact" and self.linear_mass:
            # No state-dependent dt brick. Keep the same PETSc Mat
            # rather than allocating a separate full tangent copy.
            # Do this BEFORE J/TS/KSP takes references to any Mat.
            self.tangent_mass = alias_readonly_linear_tangent(
                self.mass, self.tangent_mass
            )
        else:
            self.init_matrix(self.tangent_mass)
        self.stiffness.bindToCPU(True)
        self.init_matrix(self.stiffness)
        self.init_matrix(self.nl_stiffness)
        if self._matrix_memory_mode == "compact" and self.linear_stiffness:
            self.tangent_stiffness = alias_readonly_linear_tangent(
                self.stiffness, self.tangent_stiffness
            )
        else:
            self.init_matrix(self.tangent_stiffness)

        self.rhs = self.nl_stiffness.createVecRight()
        self.F = self.nl_stiffness.createVecRight()
        self.buffer = self.nl_stiffness.createVecRight()

        self.Scatter = PETSc.Scatter().toAll(self.F)
        self.GlobalVec = PETSc.Vec().createSeq(self.gf_model.nbdof())

    def get_cleared_TS_options(self):
        """To ensure a safe database for the PETSc TS environment"""
        
        PETSc.Options().clear()
        self.time_scheme = PETSc.Options()
        for key in self.time_scheme.getAll():
            self.time_scheme.delValue(key)

    def set_time_scheme(self, **kwargs):
        """Allows an easy setting of the PETSc TS environment

        Keyword Args:
            kwargs: PETSc TS options and SCRIMP time-stepping settings
                passed as keyword arguments (see examples).
        """

        # Options maintained by PETSc include SCRIMP-inserted defaults.
        # An entry in PETSc.Options() is NOT evidence that the caller
        # explicitly requested a particular strategy. Moreover,
        # _solve_impl() re-enters this method with zero kwargs.
        explicitly_selected_initializer = "init_step_ts_type" in kwargs
        for key, value in kwargs.items():
            self.time_scheme[key] = value

        # Linear fixed-state DAE projection is now the default
        # *candidate* for v1.3.0. Its native implicit-default
        # qualification and release-migration gates remain separate.
        # An explicitly selected legacy 'pseudo' or 'bdf' initializer
        # continues to use the original TS path.
        init_strategy = self.time_scheme.getString(
            "init_step_ts_type", "algebraic"
        ).strip().lower()
        if init_strategy not in ("algebraic", "nonlinear_affine", "pseudo", "bdf"):
            raise ValueError(
                "init_step_ts_type must be algebraic, nonlinear_affine, "
                "pseudo or bdf"
            )
        if init_strategy in ("algebraic", "nonlinear_affine") and explicitly_selected_initializer:
            # A model that deliberately disables initialization
            # must still work without specifying any strategy.
            # However, an *explicit* request for algebraic work
            # combined with init_step=false is contradictory.
            disabled = self.time_scheme.getString(
                "init_step", "true"
            ).strip().lower() in ("false", "0", "no")
            if disabled:
                raise ValueError(
                    "Explicit algebraic init_step_ts_type requires "
                    "init_step=true; omit the initializer option "
                    "when no initialization is intended"
                )

        if not self.time_scheme.hasName("ts_equation_type"):
            self.time_scheme[
                "ts_equation_type"
            ] = PETSc.TS.EquationType.DAE_IMPLICIT_INDEX2

        if not self.time_scheme.hasName("ts_type") and not self.time_scheme.hasName(
            "ts_ssp"
        ):
            self.time_scheme["ts_type"] = "bdf"
            self.time_scheme["ts_bdf_order"] = 2

        if not self.time_scheme.hasName("ksp_type"):
            self.time_scheme["ksp_type"] = "gmres"

        if not self.time_scheme.hasName("pc_type"):
            self.time_scheme["pc_type"] = "lu"

        if not self.time_scheme.hasName("pc_factor_mat_solver_type"):
            self.time_scheme["pc_factor_mat_solver_type"] = "mumps"

        # The MPI MUMPS analysis mode is the qualified v1.3 default:
        # ICNTL(28)=2 uses distributed symbolic ordering on multiple
        # ranks. Never override user-supplied options or non-MUMPS PCs.
        if (
            max_rank > 1
            and self.time_scheme.getString("pc_type", "").lower() == "lu"
            and self.time_scheme.getString(
                "pc_factor_mat_solver_type", ""
            ).lower() == "mumps"
        ):
            if not self.time_scheme.hasName("mat_mumps_icntl_28"):
                self.time_scheme["mat_mumps_icntl_28"] = 2
            if not self.time_scheme.hasName("mat_mumps_icntl_29"):
                self.time_scheme["mat_mumps_icntl_29"] = 0

        if not self.time_scheme.hasName("t_0"):
            self.time_scheme["t_0"] = 0.0

        if not self.time_scheme.hasName("t_f"):
            self.time_scheme["t_f"] = 1.0

        if not self.time_scheme.hasName("dt"):
            self.time_scheme["dt"] = 0.01

        if not self.time_scheme.hasName("dt_save"):
            self.time_scheme["dt_save"] = 0.01

        if not self.time_scheme.hasName("ts_adapt_dt_min"):
            self.time_scheme["ts_adapt_dt_min"] = 1.e-2*float(self.time_scheme["dt"])**2 / max_rank
            
        if not self.time_scheme.hasName("ts_adapt_dt_max"):
            self.time_scheme["ts_adapt_dt_max"] = float(self.time_scheme["dt_save"])

        if not self.time_scheme.hasName("ts_max_snes_failures"):
            self.time_scheme["ts_max_snes_failures"] = -1
        
        if not self.time_scheme.hasName("init_step"):
            self.time_scheme["init_step"] = True

        if not self.time_scheme.hasName("init_step_nb_iter"):
            self.time_scheme["init_step_nb_iter"] = 1

        if not self.time_scheme.hasName("init_step_ts_type"):
            # Must agree with both strategy validation above and
            # DPHS.init_step(): the previous 'pseudo' assignment
            # silently defeated the new default when no option was
            # supplied. Legacy relaxation stays explicit.
            self.time_scheme["init_step_ts_type"] = "algebraic"

        if not self.time_scheme.hasName("init_step_dt"):
            self.time_scheme["init_step_dt"] = float(self.time_scheme["dt"])**2

        self.time_scheme["isset"] = True

    def exclude_algebraic_var_from_lte(self, TS):
        """Exclude the algebraic variable from the local error troncature in the time-resolution

        Args:
            TS (PETSc TS): the PETSc TS object handling time-resolution
        """

        atol = TS.getTolerances()[1]
        atol_v = atol * np.ones((self.gf_model.nbdof(),))
        id_alg = np.array([], dtype=int)
        # We add indices of all algebraic ports
        for _, port in self.ports.items():
            if port.get_algebraic():
                I = self.gf_model.interval_of_variable(port.get_flow())
                id_alg = np.concatenate(
                    (id_alg, np.arange(I[0], I[0] + I[1], dtype=int))
                )
                # If it is not substituted, it also has the effort part
                if not port.get_substituted():
                    I = self.gf_model.interval_of_variable(port.get_effort())
                    id_alg = np.concatenate(
                        (id_alg, np.arange(I[0], I[0] + I[1], dtype=int))
                    )
            # Else on dynamical ports
            else:
                # If costate are not substituted, we also add them
                if not port.get_substituted():
                    I = self.gf_model.interval_of_variable(port.get_effort())
                    id_alg = np.concatenate(
                        (id_alg, np.arange(I[0], I[0] + I[1], dtype=int))
                    )
        id_alg.sort()
        atol_v[id_alg] = np.inf
        # Preserve the distributed layout of the TS residual vector.
        # createWithArray(global_array, comm=WORLD) would allocate N
        # local entries on EACH rank, yielding a global length P*N.
        atol_v_petsc = self.F.duplicate()
        lo, hi = atol_v_petsc.getOwnershipRange()
        atol_v_petsc.getArray()[:] = atol_v[lo:hi]
        TS.setTolerances(atol=atol_v_petsc)

    def event(self, TS, t, z, fvalue):
        """Check if the time step is not too small"""

        fvalue[0] = TS.getTimeStep() - float(self.time_scheme["ts_adapt_dt_min"])

    def postevent(self, TS, event, t, z, forward):
        """If the time step is too small, ask for the end of the simulation"""

        self.stop_TS = True
        if rank == 0:
            logging.info("\n *** AUTOMATIC STOP -- LAST ITERATION *** \n")
            logging.info(
                f"dt reached the limit {float(self.time_scheme['ts_adapt_dt_min']):8g}.\n"
            )

    def _memory_phase(self, phase, **metadata):
        """Collect a local per-rank memory sample only when enabled."""
        trace = getattr(self, "_phase_memory_trace", None)
        if trace is not None:
            trace.sample(phase, **metadata)

    def solve(self):
        """Solve the DPHS; optionally record inclusive solve timings.

        Enable SCRIMP_PROFILE=1 before creating the system, then call
        get_profile() on ALL MPI ranks to aggregate the timing report.
        """
        self.profiler.reset()
        self._profile_phase = "setup"
        # Avoid publishing a stale successful algebraic result when a
        # DPHS object is solved again with pseudo/BDF or fails early.
        self.initialization_report = None
        matrix_mode = os.environ.get(
            "SCRIMP_MATRIX_MEMORY_MODE", "legacy"
        ).strip().lower()
        if matrix_mode not in ("legacy", "compact"):
            raise ValueError(
                "SCRIMP_MATRIX_MEMORY_MODE must be legacy or compact"
            )
        self._matrix_memory_mode = matrix_mode
        # Workspaces belong to this solve, not to an old TS/KSP. Avoid
        # stale references or retaining matrices from previous solves.
        self._mass_hessian_fd_matrix = None
        self._mass_hessian_fd_probe = None
        self._mass_hessian_fd_probe_allocations = 0
        self._mass_hessian_fd_matrix_allocations = 0
        self._mass_hessian_fd_matrix_reuses = 0
        if self._explicit_rhs is not None:
            self._explicit_rhs.destroy()
            self._explicit_rhs = None
        self._explicit_rhs_ready = False
        self._explicit_rhs_time = None
        audit_dir = os.environ.get(
            "SCRIMP_MEMORY_PHASE_AUDIT_DIR", ""
        ).strip()
        self._phase_memory_trace = (
            PhaseMemoryAudit(
                audit_dir, rank=rank, ranks=max_rank,
                source_commit=os.environ.get("SCRIMP_BENCH_GIT_SHA", ""),
                mode=matrix_mode,
            )
            if audit_dir else None
        )
        self._memory_phase("solve_entry")
        backend = os.environ.get(
            "SCRIMP_TRAJECTORY_BACKEND", "memory"
        ).strip().lower()
        if backend not in ("memory", "async_binary"):
            raise ValueError(
                "SCRIMP_TRAJECTORY_BACKEND must be memory or async_binary"
            )
        self._trajectory_backend = backend
        self._trajectory_sender = None
        self._trajectory_directory = None
        if backend == "async_binary":
            from pathlib import Path
            directory = os.environ.get("SCRIMP_TRAJECTORY_DIR", "").strip()
            socket_name = os.environ.get(
                "SCRIMP_TRAJECTORY_SOCKET", ""
            ).strip()
            if not directory or not socket_name:
                raise ValueError(
                    "async_binary requires SCRIMP_TRAJECTORY_DIR "
                    "and SCRIMP_TRAJECTORY_SOCKET; start the independent "
                    "writer BEFORE mpiexec (never fork a live MPI rank)"
                )
            if len(self.solution["t"]) or len(self.solution["z"]):
                raise RuntimeError(
                    "async_binary requires an empty output trajectory"
                )
            self._trajectory_directory = Path(directory).resolve()
            if (self._trajectory_directory / "trajectory_manifest.json").exists():
                raise FileExistsError(
                    "Refusing to reuse an already committed trajectory"
                )
        # A new solve creates fresh initial/main TS/KSP/PC lifecycles:
        # never inherit numeric or symbolic snapshots from old factors.
        self._mumps_csr_patterns = {}
        self._mumps_csr_snapshots = {}
        self._mumps_reuse_stats = {}
        old_candidate = getattr(self, "_mumps_jacobian_candidate", None)
        if old_candidate is not None:
            old_candidate.destroy()
        self._mumps_jacobian_candidate = None
        self._mumps_phase_enabled = {}
        # Bricks linearity is diagnostic only; on_change now validates
        # the actual Jacobian values, including time/step dependence.
        self._all_bricks_linear = _all_bricks_are_linear(self.bricks)
        self._mumps_reset_policy = _get_mumps_reset_policy(
            all_bricks_linear=self._all_bricks_linear
        )
        _mpi_trace(
            f"MUMPS reuse policy={self._mumps_reset_policy} "
            f"(all_bricks_linear={self._all_bricks_linear})"
        )
        policies = comm.tompi4py().allgather(
            (self._mumps_reset_policy, self._matrix_memory_mode)
        )
        if any(value != policies[0] for value in policies):
            raise RuntimeError(
                "MUMPS policy and SCRIMP_MATRIX_MEMORY_MODE must resolve "
                "identically on all MPI ranks: " + repr(policies)
            )
        with self.profiler.measure("solve.total"):
            try:
                return self._solve_impl()
            finally:
                # The candidate Jacobian is only needed during callbacks.
                # Keep its allocations out of subsequent solves and do not
                # retain one extra full MPI sparse matrix after TS finishes.
                candidate = getattr(self, "_mumps_jacobian_candidate", None)
                if candidate is not None:
                    candidate.destroy()
                    self._mumps_jacobian_candidate = None
                self._mumps_csr_snapshots = {}
                self._memory_phase("solve_exit")
                if self._phase_memory_trace is not None:
                    self._phase_memory_trace.write()

    def _solve_impl(self):
        """Perform the time-resolution of the dphs thanks to PETSc TS.

        The options database is set in the `time_scheme` attribute.
        """

        assert self.time_scheme[
            "isset"
        ], "Time scheme must be defined before time-resolution"
        for name_variable in self.initial_value_set.keys():
            try:
                assert self.initial_value_set[name_variable]
            except AssertionError as err:
                logging.error(
                    f"Initial value of {name_variable} must be defined before time-resolution"
                )
                raise err
        for name_control in self.controls.keys():
            try:
                assert self.controls[name_control].get_isSet()
            except AssertionError as err:
                logging.error(
                    f"Control {name_control} must be defined before time-resolution"
                )
                raise err

        if rank == 0:
            logging.info(
                f"Simulation is starting on {comm.getSize()} processor(s) (total number of dofs: {self.gf_model.nbdof()})"
            )

        # Initiliaze the systems flags
        self.set_linear_flags()

        # Initialize time scheme options (without override)
        self.set_time_scheme()

        # Load time scheme options to PETSc
        self.time_scheme = PETSc.Options()
        self._memory_phase("before_petsc_allocate_memory")
        with self.profiler.measure("setup.petsc.allocate_memory"):
            self.allocate_memory()
        self._memory_phase("after_petsc_allocate_memory")

        # Initialization
        self.disable_all_bricks()  # enabling is done in assembly accurately
        self.gf_model.set_time(float(self.time_scheme["t_0"]))
        self.gf_model.first_iter()

        # Assembly. In compact mode a purely linear tangent aliases
        # its assembled base matrix and must not be reassembled/overwritten.
        self._memory_phase("before_getfem_assembly")
        self.assemble_mass()
        if not (self._matrix_memory_mode == "compact" and self.linear_mass):
            self.assemble_nl_mass()
        self._memory_phase("after_mass_assembly")

        self.assemble_stiffness()
        if not (self._matrix_memory_mode == "compact" and self.linear_stiffness):
            self.assemble_nl_stiffness()
        self._memory_phase("after_stiffness_assembly")

        self.assemble_rhs()
        self._memory_phase("after_rhs_assembly")

        # PETSc TS Jacobian and residual
        self.J = self.tangent_stiffness.duplicate()
        self.J.assemble()

        self.F = self.rhs.duplicate()
        self.F.assemble()

        self.buffer = self.rhs.duplicate()
        self.buffer.assemble()
        self._memory_phase("after_jacobian_and_vector_duplicates")

        self.ts_start = time.perf_counter()
        if self._checkpoint_pending is not None:
            _mpi_trace("solve: before accepted-state cold checkpoint publication")
            self._publish_pending_checkpoint_after_assembly()
            _mpi_trace("solve: verified original DAE cold checkpoint published")
        elif (
            self.time_scheme["init_step"] == "true"
        ):  # Because PETSc DB store everything as str
            _mpi_trace("solve: before initial consistency step")
            self._profile_phase = "initial"
            with self.profiler.measure("initial.consistency_total"):
                self.init_step()
            _mpi_trace("solve: after initial consistency step")
        self._memory_phase("after_initial_consistency_ts_cleanup")
        self._profile_phase = "main"

        if self._trajectory_backend == "async_binary":
            # Each producer owns its PETSc rows; no global output gather.
            local_lo, local_hi = self.F.getOwnershipRange()
            queued = int(os.environ.get(
                "SCRIMP_TRAJECTORY_MAX_QUEUED", "4"
            ))
            self._trajectory_sender = RankSender(
                os.environ["SCRIMP_TRAJECTORY_SOCKET"],
                rank=rank, ranks=max_rank,
                global_size=int(self.gf_model.nbdof()),
                lo=local_lo, hi=local_hi,
                # Match the historical monitor's real float64 archive.
                dtype=np.float64, max_queued=queued,
            )

        # TS
        TS = PETSc.TS().create(comm=comm)

        def monitor(TS, i, t, z):
            return self.monitor(
                TS,
                i,
                t,
                z,
                dt_save=float(self.time_scheme["dt_save"]),
                t_0=float(self.time_scheme["t_0"]),
            )
        
        TS.setMonitor(monitor)
        TS.setEventHandler([0], [True], self.event, self.postevent)
        TS.setEventTolerances(1e-16, vtol=[1e-19])
        TS.setIFunction(self._profile_ifunction, self.F)
        TS.setIJacobian(self._profile_ijacobian, self.J)
        if any(brick.get_explicit() for brick in self.bricks.values()):
            # PreStep runs before the Newton solve, with the last
            # accepted time/state (including the first physical step).
            TS.setPreStep(self._prestep_explicit_rhs)
        TS.setTime(float(self.time_scheme["t_0"]))
        TS.setMaxTime(float(self.time_scheme["t_f"]))
        TS.setTimeStep(float(self.time_scheme["dt"]))
        # MATCHSTEP is the default: PETSc adapts the *last* time step
        # so the integrated solution reaches t_f, rather than taking
        # a STEPOVER step past t_f and losing its final saved state.
        # A caller may explicitly restore legacy behavior using
        # ts_exact_final_time=stepover (PETSc -ts_exact_final_time).
        _set_main_ts_default_exact_final_time(TS)
        TS.setFromOptions()
        self.exclude_algebraic_var_from_lte(TS)
        InitVec = self.tangent_stiffness.createVecRight()
        InitVec.zeroEntries()
        # PETSc ownership ranges use GLOBAL indices. setValuesLocal
        # would require a local-to-global map, which is not installed.
        lo, hi = InitVec.getOwnershipRange()
        if hi > lo:
            InitVec.setValues(
                range(lo, hi),
                self.gf_model.from_variables()[lo:hi],
                addv=PETSc.InsertMode.INSERT_VALUES,
            )
        InitVec.assemble()
        
        comm.barrier()
        
        self._memory_phase("before_main_ts_solve")
        _mpi_trace("solve: before main TS.solve")
        try:
            with self.profiler.measure("main.petsc.ts_solve"):
                TS.solve(InitVec)
        except PETSc.Error:
            # Serial-only, opt-in first-step TS/SNES/KSP/PC diagnostics:
            # inspect existing solver state; never reassemble or retry.
            if (
                os.environ.get("SCRIMP_TS_FAILURE_DIAGNOSTICS") == "1"
                and max_rank == 1
            ):
                try:
                    snes = TS.getSNES()
                    ksp = snes.getKSP()
                    pc = ksp.getPC()

                    def read_state(getter):
                        try:
                            return str(getter())
                        except Exception as exc:
                            return f"unavailable:{type(exc).__name__}"

                    logging.error(
                        "MAIN TS FAILURE DIAGNOSTIC: "
                        "TS_reason=%s, TS_step=%s, TS_time=%s, "
                        "TS_dt=%s, SNES_reason=%s, SNES_iter=%s, "
                        "SNES_norm=%s, KSP_reason=%s, KSP_iter=%s, "
                        "KSP_norm=%s, PC_type=%s, PC_failed=%s, "
                        "factor_solver=%s",
                        read_state(TS.getConvergedReason),
                        read_state(TS.getStepNumber),
                        read_state(TS.getTime),
                        read_state(TS.getTimeStep),
                        read_state(snes.getConvergedReason),
                        read_state(snes.getIterationNumber),
                        read_state(snes.getFunctionNorm),
                        read_state(ksp.getConvergedReason),
                        read_state(ksp.getIterationNumber),
                        read_state(ksp.getResidualNorm),
                        read_state(pc.getType),
                        read_state(pc.getFailedReason),
                        read_state(pc.getFactorSolverType),
                    )
                except Exception as diagnostic_error:
                    logging.error(
                        "MAIN TS FAILURE diagnostic unavailable: %r",
                        diagnostic_error,
                    )
            # Requested only for a SERIAL native MUMPS failure. Read
            # KSP/PC and the actual LU factor's MUMPS status *before*
            # PETSc may destroy them. No extra assembly, no change to
            # the nonlinear PDE or the solver options. Do not attempt
            # diagnostics on only rank 0 in MPI: some getters are
            # logically collective.
            if (
                os.environ.get("SCRIMP_MUMPS_FACTOR_DIAGNOSTICS") == "1"
                and max_rank == 1
            ):
                try:
                    ksp = TS.getSNES().getKSP()
                    pc = ksp.getPC()
                    factor_details = {}
                    if pc.getType() == "lu":
                        factor = pc.getFactorMatrix()
                        factor_details = {
                            "type": factor.getType(),
                            "INFOG(1)": factor.getMumpsInfog(1),
                            "INFOG(2)": factor.getMumpsInfog(2),
                            "INFO(1)": factor.getMumpsInfo(1),
                            "INFO(2)": factor.getMumpsInfo(2),
                        } if pc.getFactorSolverType() == "mumps" else {}
                    logging.error(
                        "MAIN TS PETSc linear-solve failure: "
                        f"TS_reason={int(TS.getConvergedReason())}, "
                        f"SNES_reason={int(TS.getSNES().getConvergedReason())}, "
                        f"KSP_reason={int(ksp.getConvergedReason())}, "
                        f"PC_type={pc.getType()}, "
                        f"PC_failed={pc.getFailedReason()}, "
                        f"factor_solver={pc.getFactorSolverType()}, "
                        f"MUMPS_factor={factor_details}"
                    )
                except Exception as diagnostic_error:
                    logging.error(
                        "MAIN TS MUMPS factor diagnostic unavailable: "
                        f"{diagnostic_error!r}"
                    )
            raise
        self._memory_phase("after_main_ts_solve")
        # The stored solution samples depend on dt_save and the monitor.
        # Their last timestamp need NOT equal the actual PETSc TS end
        # time (e.g. steps at 0.0001, 0.0003 with dt_save=0.0001).
        # Keep the solver's authoritative convergence metadata BEFORE
        # TS.reset()/destroy(), for caller-side numerical acceptance.
        self.last_ts_time = float(TS.getTime())
        self.last_ts_converged_reason = int(TS.getConvergedReason())
        self.last_ts_steps = int(TS.getStepNumber())
        self.last_ts_rejections = int(TS.getStepRejections())
        self.last_ts_snes_failures = int(TS.getSNESFailures())

        # The monitor normally saves at dt_save intervals. A TS ending
        # exactly at t_f may skip the last accepted step if its advance
        # crosses a save threshold (observed for published MCSS2024:
        # accepted at 0.0001, then 0.0003, yet saved only to 0.0001).
        # Record the actual final state, without re-integrating and
        # without inventing interpolation when TS STEPOVER goes past t_f.
        t_final = float(self.time_scheme["t_f"])
        final_tol = 1e-12 * max(1.0, abs(t_final))
        if (
            self.last_ts_converged_reason > 0
            and abs(self.last_ts_time - t_final) <= final_tol
            and (
                not self.solution["t"]
                or float(self.solution["t"][-1])
                < self.last_ts_time - final_tol
            )
        ):
            self.monitor(
                TS, -1, self.last_ts_time, InitVec,
                dt_save=float(self.time_scheme["dt_save"]),
                t_0=float(self.time_scheme["t_0"]),
            )

        _mpi_trace(
            "solve: main TS completed: type={}, reason={}, "
            "t={:.9g}, dt={:.9g}, steps={}, rejected={}, "
            "SNES_failures={}".format(
                TS.getType(), int(TS.getConvergedReason()),
                float(TS.getTime()), float(TS.getTimeStep()),
                TS.getStepNumber(), TS.getStepRejections(),
                TS.getSNESFailures(),
            )
        )
        del InitVec

        if rank == 0:
            logging.info(f"Elapsed time: {time.perf_counter() - self.ts_start:1.4g}s")
            logging.info(
                f"Steps: {TS.getStepNumber()} ({TS.getStepRejections()} rejected, {TS.getSNESFailures()} Nonlinear solver failures)"
            )
            logging.info(
                f"Nonlinear iterations: {TS.getSNESIterations()}, Linear iterations: {TS.getKSPIterations()}"
            )
        
        with self.profiler.measure("main.petsc.cleanup"):
            TS.reset()
            TS.destroy()
            if self._explicit_rhs is not None:
                self._explicit_rhs.destroy()
                self._explicit_rhs = None
            self._explicit_rhs_ready = False
            self._explicit_rhs_time = None
            gc.collect()
        self._memory_phase("after_main_ts_cleanup")

        # Preserve the real end time and publish the trajectory only
        # after all local transport threads have drained. A committed
        # manifest is mandatory; no incomplete output can look solved.
        if self._trajectory_backend == "async_binary" and (
            self.last_ts_converged_reason <= 0
        ):
            raise RuntimeError(
                "Refusing to commit an unconverged PETSc trajectory"
            )
        if self._trajectory_backend == "async_binary":
            self._trajectory_sender.close()
            self._trajectory_sender = None
            manifest = wait_for_manifest(
                self._trajectory_directory / "trajectory_manifest.json",
                timeout=float(os.environ.get(
                    "SCRIMP_TRAJECTORY_COMMIT_TIMEOUT_S", "30"
                )),
            )
            saved_t, indices = (
                _select_saved_time_indices_through_final_time(
                    self.solution["t"], float(self.time_scheme["t_f"])
                )
            )
            archived = LazyTrajectory(manifest, indices=indices)
            raw_t = np.asarray(self.solution["t"], dtype=float)
            if (len(archived) != len(saved_t)
                    or not np.array_equal(archived.times, raw_t[indices])):
                raise RuntimeError(
                    "Binary archive metadata differs from PETSc monitor"
                )
            self.solution["t"] = saved_t
            self.solution["z"] = archived
        else:
            self.solution["t"], self.solution["z"] = (
                _select_saved_solution_through_final_time(
                    self.solution["t"], self.solution["z"],
                    float(self.time_scheme["t_f"]),
                )
            )

        self.enable_all_bricks()
        self.solve_done = True

    def _sync_accepted_pseudo_initial_state(self, init_vec):
        """Publish the accepted PETSc pseudo-TS state to GetFEM on all ranks.

        A DAE pseudo-initialization can change both algebraic and
        differential unknowns. The following physical TS reconstructs
        its initial Vec from gf_model.from_variables(), so merely
        solving the pseudo-TS does NOT propagate its accepted state.

        In v1.2.0 the monitor performed this synchronization only as
        a side effect of linear_stiffness=False; prescribed source
        bricks happened to set that flag. The correct classification
        of RHS-only sources must not suppress consistent initialization.
        """
        with self.profiler.measure("initial.petsc.accepted_state_sync"):
            publish_accepted_pseudo_state(
                self.Scatter[0], init_vec, self.GlobalVec, self.gf_model
            )

    def init_step(self):
        """Establish a DAE-consistent physical initial state.

        The v1.3.0 default candidate `algebraic` projects only
        algebraic variables while preserving prescribed differential
        initial coefficients. It is linear-only; unsupported or
        inconsistent initial data raise explicit errors.
        Historical `pseudo` remains available by explicit request;
        it integrates a tiny pseudo-time PETSc TS and can move the
        differential coefficients. The alternative historical
        `bdf` also uses a short time integration, optionally from a
        negative initial time with `init_step_bdf_order=2` and
        `init_step_from_negative_time=True`. PETSc BDF2 bootstraps
        its history internally with a half-step (PETSc 3.21).

        `init_step_ts_type=algebraic` directly projects only
        algebraic variables with a distributed PETSc LSQR solve, and
        performs no pseudo-time integration. Explicit historical
        relaxation options remain available. Changing the default
        is subject to an independent native and migration release gate;
        the main physical TS remains unchanged.

        See docs/V130_DAE_ALGEBRAIC_INITIALIZER_DESIGN.md.
        """

        init_type = self.time_scheme.getString(
            "init_step_ts_type", "algebraic"
        ).strip().lower()
        if init_type == "algebraic":
            with self.profiler.measure("initial.petsc.algebraic_projection"):
                self.initialization_report = initialize_linear_algebraic(self)
            if rank == 0:
                result = self.initialization_report
                logging.info(
                    "Algebraic DAE initial projection: "
                    "true ||F_a||=%.9g, feasibility tol=%.9g, "
                    "accuracy goal=%.9g, accuracy_met=%s, passes=%d",
                    result["native_final_residual_l2"],
                    result["feasibility_tolerance"],
                    result["accuracy_target"],
                    result["accuracy_target_met"],
                    len(result["passes"]),
                )
            return
        if init_type == "nonlinear_affine":
            # Opt-in ONLY: full model may be nonlinear in h,p and in
            # the mass. Algebraic block must be regular and affine in
            # algebraic coefficients at original fixed differential q0.
            # Singular high-index boundary multipliers are rejected.
            with self.profiler.measure("initial.petsc.nonlinear_affine_projection"):
                self.initialization_report = (
                    initialize_nonlinear_affine_algebraic(self)
                )
            if rank == 0:
                logging.info(
                    "Nonlinear-state affine-algebraic initial projection: "
                    "native ||F_a||=%.9g, target=%.9g",
                    self.initialization_report["native_final_algebraic_residual_l2"],
                    self.initialization_report["accuracy_target"],
                )
            return
        if init_type not in ("pseudo", "bdf"):
            raise ValueError(
                "init_step_ts_type must be algebraic, nonlinear_affine, "
                "pseudo or bdf"
            )
        init_dt = float(self.time_scheme["init_step_dt"])
        if not np.isfinite(init_dt) or init_dt <= 0:
            raise ValueError("init_step_dt must be finite and positive")
        negative_bdf_start = (
            init_type == "bdf"
            and self.time_scheme.getString(
                "init_step_from_negative_time", "false"
            ).lower() == "true"
        )
        t0 = float(self.time_scheme["t_0"])

        # TS
        TS = PETSc.TS().create(comm=comm)

        def monitor(TS, i, t, z):
            return self.monitor(
                TS,
                i,
                t,
                z,
                dt_save=1.0,
                t_0=float(self.time_scheme["t_0"]),
                initial_step=True,
            )

        TS.setMonitor(monitor)
        TS.setIFunction(self._profile_ifunction, self.F)
        TS.setIJacobian(self._profile_ijacobian, self.J)
        if any(brick.get_explicit() for brick in self.bricks.values()):
            # PreStep runs before the Newton solve, with the last
            # accepted time/state (including the first physical step).
            TS.setPreStep(self._prestep_explicit_rhs)
        TS.setTime(t0 - init_dt if negative_bdf_start else t0)
        TS.setMaxSteps(
            1 if negative_bdf_start
            else int(self.time_scheme["init_step_nb_iter"])
        )
        if negative_bdf_start:
            # A single full step ends exactly at the physical t0. PETSc
            # initializes the order-2 BDF history at t0-dt/2.
            TS.setMaxTime(t0)
        TS.setExactFinalTime(PETSc.TS.ExactFinalTime.MATCHSTEP)

        # Save options for later use and avoid overidden them with the initial step
        saved_options = dict()
        if self.time_scheme.hasName("ts_type"):
            saved_options["ts_type"] = self.time_scheme["ts_type"]
            self.time_scheme.delValue("ts_type")
        if self.time_scheme.hasName("ts_bdf_order"):
            saved_options["ts_bdf_order"] = self.time_scheme["ts_bdf_order"]
            self.time_scheme.delValue("ts_bdf_order")
        if self.time_scheme.hasName("ts_adapt_type"):
            saved_options["ts_adapt_type"] = self.time_scheme["ts_adapt_type"]
            if negative_bdf_start:
                self.time_scheme.delValue("ts_adapt_type")
        if self.time_scheme.hasName("ts_ssp"):
            saved_options["ts_ssp"] = self.time_scheme["ts_ssp"]
            self.time_scheme.delValue("ts_ssp")
        if self.time_scheme.hasName("dt"):
            saved_options["dt"] = float(self.time_scheme["dt"])
        if self.time_scheme.hasName("ts_adapt_dt_min"):
            saved_options["ts_adapt_dt_min"] = float(
                self.time_scheme["ts_adapt_dt_min"]
            )
            # Only override the pseudo-TS minimum step temporarily.
            # Deleting dt_max here was a typo: it could accidentally
            # alter the initial TS adaptation configuration.
            self.time_scheme.delValue("ts_adapt_dt_min")
        if self.time_scheme.hasName("ts_adapt_dt_max"):
            saved_options["ts_adapt_dt_max"] = float(
                self.time_scheme["ts_adapt_dt_max"]
            )
            self.time_scheme.delValue("ts_adapt_dt_max")
            
        self.time_scheme["ts_type"] = init_type
        if init_type == "bdf":
            self.time_scheme["ts_bdf_order"] = int(
                self.time_scheme.getString("init_step_bdf_order", "1")
            )
        if negative_bdf_start:
            # A fixed step is essential: otherwise TSAdapt may move
            # the endpoint away from the intended physical t0.
            self.time_scheme["ts_adapt_type"] = "none"
        self.time_scheme["ts_adapt_dt_min"] = 1e-24

        TS.setTimeStep(init_dt)
        
        _mpi_trace("init_step: before TS.setFromOptions")
        TS.setFromOptions()

        # Diagnostic-only alternative to the global LU/MUMPS factorization.
        # The default solver remains unchanged. The native MPI crash observed
        # in issue #6 occurs in MatDestroy_MUMPS -> PCReset_LU -> TSReset.
        # Using BJACOBI in the initial pseudo-TS alone helps distinguish a
        # MUMPS-specific factor lifecycle defect from an earlier heap write.
        init_pc_type = os.environ.get("SCRIMP_INIT_PC_TYPE", "").strip()
        init_solver = os.environ.get(
            "SCRIMP_INIT_FACTOR_SOLVER_TYPE", ""
        ).strip().lower()
        if init_pc_type or init_solver:
            init_ksp = TS.getSNES().getKSP()
            init_pc = init_ksp.getPC()
            if init_pc_type:
                if init_pc_type not in ("bjacobi", "none", "lu"):
                    raise ValueError(
                        "SCRIMP_INIT_PC_TYPE must be bjacobi, none or lu"
                    )
                init_pc.setType(init_pc_type)
                if init_pc_type in ("bjacobi", "none"):
                    # KSPPREONLY would apply just one preconditioner
                    # action, not solve the system. Use GMRES.
                    init_ksp.setType("gmres")
                    init_ksp.setTolerances(
                        rtol=1e-8, atol=1e-12, max_it=500
                    )
            if init_solver:
                # Diagnostically mix TWO LU packages on the SAME
                # published FE model: e.g. SuperLU only for the
                # consistent-state pseudo-TS; MUMPS for main BDF2.
                # This does NOT change the global PETSc options or
                # the main TS's factor package.
                if init_solver not in ("superlu", "mumps"):
                    raise ValueError(
                        "SCRIMP_INIT_FACTOR_SOLVER_TYPE must be "
                        "superlu or mumps"
                    )
                if init_pc.getType() != "lu":
                    raise ValueError(
                        "SCRIMP_INIT_FACTOR_SOLVER_TYPE requires "
                        "pc_type=lu during initialization"
                    )
                init_pc.setFactorSolverType(init_solver)
            if rank == 0:
                logging.info(
                    "Initial consistency solver override: "
                    f"KSP={init_ksp.getType()}, PC={init_pc.getType()}, "
                    f"LU={init_solver or 'inherited'} (initial TS only)"
                )

        _mpi_trace("init_step: before algebraic tolerance setup")
        self.exclude_algebraic_var_from_lte(TS)
        _mpi_trace("init_step: TS setup completed")

        if rank == 0:
            logging.info(
                "Initial consistency: TS={} (BDF order={}), steps={}, "
                "start={:.12g}, target={:.12g}, dt={:.12g}, "
                "adapt={}".format(
                    init_type,
                    self.time_scheme.getString("ts_bdf_order", "n/a"),
                    1 if negative_bdf_start else int(
                        self.time_scheme["init_step_nb_iter"]
                    ),
                    t0 - init_dt if negative_bdf_start else t0,
                    t0 if negative_bdf_start else float("nan"),
                    init_dt,
                    self.time_scheme.getString("ts_adapt_type", "default"),
                )
            )

        InitVec = self.tangent_stiffness.createVecRight()
        InitVec.zeroEntries()
        # PETSc ownership ranges use GLOBAL indices. setValuesLocal
        # would require a local-to-global map, which is not installed.
        lo, hi = InitVec.getOwnershipRange()
        if hi > lo:
            InitVec.setValues(
                range(lo, hi),
                self.gf_model.from_variables()[lo:hi],
                addv=PETSc.InsertMode.INSERT_VALUES,
            )
        InitVec.assemble()
        initial_state_norm = float(InitVec.norm())
        # Optional diagnostics: distinguish an algebraic consistency solve
        # from pseudo-time relaxation that changes the physical h and p.
        # Compare the actual PETSc final state, not the last Newton trial
        # held by GetFEM. Copies are only allocated for this opt-in audit.
        audit_dynamic = os.environ.get(
            "SCRIMP_INIT_DYNAMIC_STATE_AUDIT", "0"
        ) == "1"
        if audit_dynamic:
            initial_dynamics = {}
            initial_all = np.asarray(
                self.gf_model.from_variables()
            ).real.copy()
            for name in self.states:
                start, width = self.gf_model.interval_of_variable(name)
                start, width = int(start), int(width)
                initial_dynamics[name] = initial_all[
                    start:start + width
                ].copy()
        comm.barrier()

        _mpi_trace("init_step: before consistency TS.solve")
        try:
            with self.profiler.measure("initial.petsc.ts_solve"):
                TS.solve(InitVec)
        except PETSc.Error:
            # A TS_DIVERGED_STEP_REJECTED is not an MUMPS error.
            # Expose the nested SNES/KSP reasons before PETSc's
            # traceback; a failure of the nonlinear convergence
            # test differs fundamentally from failed LU factorization.
            if rank == 0:
                snes = TS.getSNES()
                ksp = snes.getKSP()
                logging.error(
                    "Initial consistency FAILED: "
                    f"TS_reason={int(TS.getConvergedReason())}, "
                    f"steps={TS.getStepNumber()}, "
                    f"rejected={TS.getStepRejections()}, "
                    f"SNES_failures={TS.getSNESFailures()}, "
                    f"SNES_reason={int(snes.getConvergedReason())}, "
                    f"SNES_iterations={snes.getIterationNumber()}, "
                    f"KSP_reason={int(ksp.getConvergedReason())}, "
                    f"KSP_iterations={ksp.getIterationNumber()}"
                )
            raise
        final_state_norm = float(InitVec.norm())

        # The accepted PETSc pseudo-TS solution is the authoritative
        # consistent initial state for the main TS. Previously, its
        # propagation to GetFEM depended accidentally on whether a
        # source brick marked linear_stiffness=False: the monitor
        # then called to_variables(), but correctly classified
        # prescribed RHS-only sources now leave linear_stiffness=True.
        # Synchronize explicitly for ALL models (including linear,
        # nonlinear, MPI and lumped PDE/ODE), independently of those
        # Jacobian flags. Do this before TS.reset() and before the main
        # InitVec is reconstructed from gf_model.from_variables().
        # This collective scatter populates GlobalVec on every rank.
        if int(TS.getConvergedReason()) > 0 and TS.getStepNumber() > 0:
            self._sync_accepted_pseudo_initial_state(InitVec)
        if audit_dynamic:
            # Compare the exact accepted pseudo-TS state, now also
            # synchronized to GetFEM for the following main solve.
            if int(TS.getConvergedReason()) <= 0 or TS.getStepNumber() <= 0:
                self.Scatter[0].scatter(InitVec, self.GlobalVec)
            final_all = np.asarray(
                self.GlobalVec.getArray(readonly=True)
            ).real
            if rank == 0:
                for name, initial_values in initial_dynamics.items():
                    start, width = self.gf_model.interval_of_variable(name)
                    final_values = final_all[
                        int(start):int(start + width)
                    ]
                    delta = float(np.linalg.norm(
                        final_values - initial_values
                    ))
                    baseline = float(np.linalg.norm(initial_values))
                    print(
                        f"[SCRIMP initial dynamic audit] {name}: "
                        f"initial_norm={baseline:.9g}, "
                        f"final_norm={np.linalg.norm(final_values):.9g}, "
                        f"absolute_change={delta:.9g}, "
                        f"relative_change={delta / max(baseline, 1e-14):.9g}",
                        file=sys.stderr, flush=True,
                    )
        _mpi_trace(
            "init_step: consistency TS completed: reason={}, "
            "t={:.12g}, steps={}, rejected={}, SNES_failures={}, "
            "initial_norm={:.9g}, final_norm={:.9g}".format(
                int(TS.getConvergedReason()), float(TS.getTime()),
                TS.getStepNumber(), TS.getStepRejections(),
                TS.getSNESFailures(),
                initial_state_norm, final_state_norm,
            )
        )
        if negative_bdf_start:
            if int(TS.getConvergedReason()) <= 0 or (
                abs(float(TS.getTime()) - t0)
                > 1e-12 * max(1.0, abs(t0))
            ):
                raise RuntimeError(
                    "Initial BDF2 consistency solve failed to reach "
                    f"t0={t0:.12g}: reached {TS.getTime()}, "
                    f"reason={TS.getConvergedReason()}"
                )
        # The norm guard is opt-in through init_step_max_norm_amplification.
        # Apply it to *any* initial TS type when set: previously it only
        # protected negative-start BDF2 and could silently accept the
        # disastrous ~1e24 pseudo-TS result in the SWE diagnostic.
        # Passing this guard does NOT prove DAE consistency or that
        # the original dynamical initial states have been preserved.
        if self.time_scheme.hasName("init_step_max_norm_amplification"):
            norm_limit = float(self.time_scheme[
                "init_step_max_norm_amplification"
            ])
            if np.isnan(norm_limit) or norm_limit <= 0:
                raise ValueError(
                    "init_step_max_norm_amplification must be positive"
                )
            if int(TS.getConvergedReason()) <= 0:
                raise RuntimeError(
                    f"Initial {init_type} TS failed: "
                    f"reason={int(TS.getConvergedReason())}, "
                    f"steps={TS.getStepNumber()}"
                )
            if not np.isfinite(final_state_norm) or (
                final_state_norm
                > norm_limit * max(initial_state_norm, 1.0)
            ):
                label = "BDF2" if negative_bdf_start else init_type.upper()
                raise RuntimeError(
                    f"Initial {label} state diverged: "
                    f"||z_initial||={initial_state_norm:.9g}, "
                    f"||z_final||={final_state_norm:.9g}, "
                    f"maximum factor={norm_limit:.9g}"
                )
        del InitVec

        if rank == 0:
            logging.info(
                f"Initialisation done in {time.perf_counter() - self.ts_start:1.4g}s"
            )

        _mpi_trace("init_step: before TS.reset (PC factor cleanup)")
        with self.profiler.measure("initial.petsc.cleanup"):
            TS.reset()
            _mpi_trace("init_step: after TS.reset")
            TS.destroy()
            _mpi_trace("init_step: after TS.destroy")
            gc.collect()

        # Delete temporary options, including BDF2/adapt overrides
        # which must not leak into the main time integrator.
        for key in (
            "ts_type", "ts_bdf_order", "ts_adapt_type",
            "ts_adapt_dt_min",
        ):
            if self.time_scheme.hasName(key):
                self.time_scheme.delValue(key)

        # Restore the main TS options exactly as they were.
        for key in saved_options.keys():
            self.time_scheme[key] = saved_options[key]

    def monitor(self, TS, i, t, z, dt_save=1.0, t_0=0.0, initial_step=False):
        """Monitor to use during iterations of time-integration at each successful time step

        Args:
            TS (PETSc TS): the PETSc TS object handling time-resolution
            i (int): the iteration in the time-resolution
            t (float): time parameter
            z (PETSc Vec): the state
            dt_save (float): save the solution each dt_save s (different from the time-step `dt` used for resolution)
            t_0 (float): the initial time
            initial_step (bool): `True` if this is the initial consistency step (default=`False`)
        """

        # Step 0 precedes the main factorization; step 1 follows it.
        # Sampling is local-only; it never introduces MPI collectives.
        if not initial_step and i in (0, 1):
            self._memory_phase(
                "main_monitor_step_" + str(i),
                step=int(i), time=float(t),
            )

        # Do not gather full state on all ranks merely for archival.
        # The nonlinear GetFEM model STILL needs its global transient
        # state; PETSc's Newton trial-state path is unchanged.
        needs_global_state = (
            initial_step
            or getattr(self, "_trajectory_backend", "memory") != "async_binary"
            or not self.linear_mass or not self.linear_stiffness
        )
        if needs_global_state:
            self.Scatter[0].scatter(z, self.GlobalVec)

        if not initial_step:
            if self.solution["t"]:
                next_saved_at_t = self.solution["t"][-1] + dt_save
            else:
                next_saved_at_t = t_0
                
            spent_time = time.perf_counter()-self.ts_start
            current_dt = float(TS.getTimeStep())
            remaining_time = ((float(self.time_scheme["t_f"])+current_dt)/(t+current_dt)-1) * spent_time
            progress = int(100*t/float(self.time_scheme["t_f"]))
            # Store the exact final state of coupled PDE/ODE runs even if
            # the next dt_save threshold would otherwise skip it.
            final_lumped_step = (
                bool(self.lumped_variables)
                and t >= float(self.time_scheme["t_f"]) - 1.e-12 * max(
                    1.0, abs(float(self.time_scheme["t_f"]))
                )
            )
            if (
                (i == 0) or (next_saved_at_t - t <= 0)
                or (i == -1) or self.stop_TS or final_lumped_step
            ):
                if rank == 0:
                    sys.stdout.write(
                        f"\ri={i:8d} t={t:8g} * ({int(spent_time)}s)   dt={current_dt:8g}   ETE={int(remaining_time)}s   ({progress}%)        \n"
                    )
                    sys.stdout.flush()
                if getattr(self, "_trajectory_backend", "memory") == "async_binary":
                    # Short rank-local copy only. Disk/socket waits happen
                    # in the I/O thread, never in this TS callback.
                    self._trajectory_sender.enqueue(
                        i, t, z.getArray(readonly=True).real
                    )
                    self.solution["t"].append(t)
                else:
                    self.solution["t"].append(t)
                    self.solution["z"].append(
                        np.array(
                            self.GlobalVec.getArray(readonly=True).real,
                            dtype=np.float64, copy=True,
                        )
                    )
            else:
                if rank == 0:
                    sys.stdout.write(
                        f"\ri={i:8d} t={t:8g}   ({int(spent_time)}s)   dt={current_dt:8g}   ETE={int(remaining_time)}s   ({progress}%)          "
                    )
                    sys.stdout.flush()

        comm.barrier()
        self.gf_model.set_time(t)  # Update the time t in the getfem `Model`
        if (
            not self.linear_mass or not self.linear_stiffness
        ):  # Needed for non-linear systems only
            # Update the state in the getfem `Model`
            self.gf_model.to_variables(
                np.array(
                    self.GlobalVec.getArray(readonly=True).real,
                    dtype=np.float64, copy=True,
                )
            )
        # Says to getfem that we go to the next iteration, not sure if needed
        self.gf_model.next_iter()

        # Do NOT force PETSc.garbage_cleanup() inside a TS callback.
        # The former call triggered PetscGarbageCleanup -> MatDestroy_MPIAIJ
        # -> free(): invalid pointer in a 2-rank pseudo-TS run. PETSc/TS
        # owns live matrix references here; leave object cleanup outside
        # the solver and avoid transient Mat allocations in Newton calls.

    def compute_Hamiltonian(self):
        """Compute each `term` constituting the Hamiltonian

        Wrapper to the function `compute` of Hamiltonian class
        """

        try:
            assert self.solve_done
        except AssertionError as err:
            logging.error(
                "System has not been solved yet, Hamiltonian can not be computed"
            )
            raise err

        self.hamiltonian.compute(self.solution, self.gf_model, self.domain)
    
    def get_Hamiltonian(self):
        """This function returns the Hamiltonian in function of time
        
        Returns:
            numpy array: time array of the values of the Hamiltonian.
        """
        
        if not self.hamiltonian.get_is_computed():
            self.compute_Hamiltonian()

        t = np.array(self.solution["t"])
        HamTot = np.zeros(t.size)
        for _, term in enumerate(self.hamiltonian):
            values = np.array(term.get_values())
            HamTot += values
            
        return HamTot

    def plot_Hamiltonian(self, with_powers=True, save_figure=False, filename="Hamiltonian.png"):
        """Plot each term constituting the Hamiltonian and the Hamiltonian

        May include the `power terms`, i.e. the sum over [t_0, t_f] of the flow/effort product of algebraic ports.

        In MPI runs this method must be called on every rank, since computing
        Hamiltonian terms and port powers requires collective GetFEM assembly.
        Only rank 0 creates and displays the plot.

        Args:
            - with_powers (bool): if `True` (default), the plot will also contains the power of each algebraic ports
            - save_figure (bool): if 'True' (defaults: False), save the plot
            - filename (str): the name of the file where the plot is saved (defaults: `Hamiltonian.png`)
        """

        if not self.hamiltonian.get_is_computed():
            self.compute_Hamiltonian()

        t = np.array(self.solution["t"])
        HamTot = np.zeros(t.size)
        ax = None
        if rank == 0:
            fig = plt.figure(figsize=[8, 5])
            ax = fig.add_subplot(111)

        for _, term in enumerate(self.hamiltonian):
            values = np.array(term.get_values())
            HamTot += values
            if rank == 0:
                ax.plot(t, values, label=term.get_description())

        if len(self.hamiltonian) > 1 and rank == 0:
            ax.plot(t, HamTot, label=self.hamiltonian.get_name())

        # GetFEM power assembly can involve MPI collectives: every rank must
        # enter plot_powers(), even though only rank 0 draws the result.
        if with_powers:
            self.plot_powers(ax, HamTot=HamTot)

        if rank == 0:
            ax.legend()
            ax.grid(axis="both")
            ax.set_xlabel("time t")
            ax.set_ylabel("Hamiltonian terms")
            ax.set_title("Evolution of Hamiltonian terms")
            if save_figure:
                plt.savefig(os.path.join(outputs_path, "png", filename), dpi=300)
            plt.show()

    def compute_powers(self):
        """Compute raw signed algebraic-port powers with typed sampling.

        For CN, preserve the already validated step-associated GetFEM
        powers exactly. For BDF/ARKIMEX/other schemes, compute generic
        instantaneous powers on saved accepted states; PETSc's hidden
        internal stages are NOT reconstructed from two endpoints.
        The public Port.get_power_sampling() documents the choice.
        """

        try:
            assert self.solve_done
        except AssertionError as err:
            logging.error("System has not been solved yet, powers can not be computed")
            raise err

        if rank == 0:
            logging.info(
                "Start computing the powers (substituted ports are not automated)"
            )
        
        start = time.perf_counter()
        ts_type = self.time_scheme.getString("ts_type", "bdf").lower()
        sampling = "cn_step" if ts_type == "cn" else "instantaneous"
        self.powers_computed = False
        for _, port in self.ports.items():
            port.compute(
                self.solution, self.gf_model, self.domain, sampling=sampling
            )

        self.powers_computed = True

        if rank == 0:
            logging.info(
                f"Powers have been computed in {time.perf_counter() - start} s"
            )

    def plot_powers(self, ax=None, HamTot=None):
        """Plot each power associated to each algebraic port

        The numerical quadrature follows each Port's declared
        get_power_sampling() convention: 'cn_step' integrates original
        step powers P[k>=1] as sum(dt[k]*P[k]), never reusing
        instantaneous P[0]; 'instantaneous' integrates actual nodal
        P(t[k],z[k]) with a trapezoidal *approximation*. An unknown
        convention fails explicitly, preventing a silent wrong plot.

        CN is already native-qualified for the tested mixed-boundary
        wave; exact BDF/ARKIMEX stage quadrature is NOT claimed
        without the underlying PETSc internal-stage data. Port signs
        are the original raw flow-effort orientation.
        If HamTot is provided, the plotted Balance retains the historic
        sign convention; its constancy is model-dependent.

        ax (Matplotlib axis): the gca of matplotlib when this function
            is called from plot_Hamiltonian().
        HamTot (Numpy array): the values of the Hamiltonian over time
        """

        if not self.powers_computed:
            self.compute_powers()

        need_to_set_figure = False

        t = np.array(self.solution["t"])
        if ax is None and rank == 0:
            fig = plt.figure(figsize=[8, 5])
            ax = fig.add_subplot(111)
            need_to_set_figure = True

        if HamTot is not None:
            SP_balance = HamTot.copy()

        for _, port in self.ports.items():
            if port.get_algebraic() and not port.get_substituted():
                power = port.get_power()
                sampling = port.get_power_sampling()
                if sampling == "cn_step":
                    convention = "cn_step"
                    display_convention = "CN step integral"
                elif sampling == "instantaneous":
                    convention = "instantaneous_trapezoid"
                    display_convention = "nodal trapezoid (approximate)"
                else:
                    raise ValueError(
                        "Cannot integrate port " + port.get_name()
                        + ": unknown power sampling convention"
                    )
                values = integrate_port_power(
                    t, power, convention=convention
                )
                if rank == 0:
                    if port.get_dissipative():
                        ax.plot(t, -values, label=port.get_name()
                                + " (" + display_convention + ")")
                    else:
                        ax.plot(t, values, label=port.get_name()
                                + " (" + display_convention + ")")
                if HamTot is not None:
                    if port.get_dissipative():
                        SP_balance += values
                    else:
                        SP_balance -= values

        if HamTot is not None:
            # Check if the balance makes sense: should not have a port which is both algebraic and substituted
            check_makes_sense = True
            for _, port in self.ports.items():
            	if port.get_algebraic() and port.get_substituted():
                    check_makes_sense = False
            if check_makes_sense and rank == 0:
                ax.plot(t, SP_balance, "--",
                        label="Balance (legacy port signs; model-dependent)")

        if need_to_set_figure and rank == 0:
            ax.legend()
            ax.grid(axis="both")
            ax.set_xlabel("time t")
            ax.set_ylabel("Integrated signed port powers")
            ax.set_title("Cumulative power terms (sampling explicit)")
            plt.show()

    def spy_Dirac(self, t=None, state=None):
        """
        !TO DO: improve a lot with position of bricks and no constitutive relations!!!
        """

        # TODO: improve a lot!!!

        if not t == None:
            self.gf_model.set_time(t)
        else:
            self.gf_model.set_time(0.0)

        if not state == None:
            self.gf_model.to_variables(state)
        else:
            self.gf_model.to_variables(np.zeros((1, self.gf_model.nbdof())))

        size = self.gf_model.nbdof()
        M = 0
        J = 0

        # Avoid missing a brick!
        self.enable_all_bricks()

        for _, brick in self.bricks.items():
            if brick.get_position() == "flow":
                for region in brick.get_regions():
                    matrix = gf.asm_generic(
                        self.domain._int_method[brick.get_mesh_id()],
                        2,
                        brick.get_form(),
                        region,
                        self.gf_model,
                    )
                    M += extract_gmm_to_scipy([0, size], [0, size], matrix)
            if brick.get_position() == "effort":
                for region in brick.get_regions():
                    matrix = gf.asm_generic(
                        self.domain._int_method[brick.get_mesh_id()],
                        2,
                        brick.get_form(),
                        region,
                        self.gf_model,
                    )
                    J += extract_gmm_to_scipy([0, size], [0, size], matrix)

        pl.spy(M, markersize=0.2)
        pl.show()
        pl.spy(J, markersize=0.2)
        pl.show()

    def export_matrices(self, t=None, state=None, path=None, to="matlab"):
        """
        TODO:
        """

        # TODO

        if path == None:
            path = outputs_path

    def export_to_pv(self, name_variable, path=None, t="All"):
        """Export the solution to .vtu file(s) (for ParaView), with associated .pvd if t='All'

        name_variable (str): the variable to export
        path (str): the path for the output file (default: in the `outputs/pv` folder next to your .py file)
        t (Numpy array): the time values of extraction (default: `All` the stored times), `All`, `Init`, `Final`
        """

        try:
            assert self.solve_done
        except AssertionError as err:
            logging.error("The system has not been solved yet")
            raise err

        if path is not None:
            path = os.path.join(path, "pv", name_variable)
        else:
            path = os.path.join(outputs_path, "pv", name_variable)

        if rank == 0:
            # Old pv must be removed first
            if os.path.exists(path):
                import shutil
                shutil.rmtree(path, ignore_errors=True)
            os.makedirs(path)
        comm.barrier()

        # Get the dofs of name_variable
        I = self.gf_model.interval_of_variable(name_variable)
        J = range(I[0], I[0] + I[1])

        if t == "All":
            if rank == 0:
                logging.info(
                    f"Export all time values of {name_variable} is starting..."
                )
            with open(os.path.join(path, name_variable + ".pvd"), "w") as pvd_file:
                pvd_file.write('<?xml version="1.0"?>\n')
                pvd_file.write('<VTKFile type="Collection" version="0.1">\n')
                pvd_file.write("  <Collection>\n")
                for k, _ in enumerate(self.solution["t"]):
                    z = self.solution["z"][k][J]
                    FEM = self.gf_model.mesh_fem_of_variable(name_variable)
                    FEM.export_to_vtu(
                        os.path.join(path, name_variable + f"_{k:05d}.vtu"),
                        FEM,
                        z,
                        name_variable,
                    )
                    pvd_file.write(
                        '    <DataSet timestep="'
                        + str(self.solution["t"][k])
                        + '" file="'
                        + name_variable
                        + f"_{k:05d}.vtu"
                        + '" />\n'
                    )
                    if rank == 0:
                        sys.stdout.write("\rExport time %f" % self.solution["t"][k])
                        sys.stdout.flush()
                pvd_file.write("  </Collection>\n")
                pvd_file.write("</VTKFile>")
            pvd_file.close()
            if rank == 0:
                logging.info("\rExport is done              ")
        elif t == "Init":
            if rank == 0:
                logging.info(f"Export initial value of {name_variable}")
            z = self.solution["z"][0][J]
            FEM = self.gf_model.mesh_fem_of_variable(name_variable)
            FEM.export_to_vtu(
                os.path.join(path, name_variable + "_init.vtu"),
                FEM,
                z,
                name_variable,
            )
        elif t == "Final":
            if rank == 0:
                logging.info(f"Export final value of {name_variable}")
            z = self.solution["z"][-1][J]
            FEM = self.gf_model.mesh_fem_of_variable(name_variable)
            FEM.export_to_vtu(
                os.path.join(path, name_variable + "_final.vtu"),
                FEM,
                z,
                name_variable,
            )
        else:
            logging.error(f"t must be `All`, `Init` or `Final`. Unknown value: {t}")
            raise ValueError

    def get_solution(self, name_variable) -> list:
        """This functions is useful (especially in 1D) to handle post-processing
        by extracting the variable of interest from the PETSc Vec self.solution["z"]

        Args:
            name_variable (str): the name of the variable to extract

        Returns:
            list(numpy array): a list of numpy array (dofs of name_variable) at each time step (according to self.solution["t"])
        """

        try:
            assert self.solve_done
        except AssertionError as err:
            logging.error("The system has not been solved yet")
            raise err

        dofs = self.gf_model.interval_of_variable(name_variable)
        variable = []
        for k, _ in enumerate(self.solution["t"]):
            variable.append(
                np.array(self.solution["z"][k][dofs[0] : dofs[0] + dofs[1]])
            )

        return variable

    def get_quantity(self, expression, region=-1, order=0, CN=False, mesh_id=0) -> list:
        """This functions computes the integral over region of the expression
        at each time step

        Args:
            - expression (str): the GFWL expression to compute
            - region (int, optional): the id of the region (defaults to -1)
            - order (int, optional): the order of the quantity to be computed (0: scalar, 1: vector, 2: tensor) (defaults to 0)
            - CN (bool, optional): take half the sum 0.5*(z^{n+1}+z^n) instead of z^n in the computation (defaults to False)
            - mesh_id (int, optional): the id of the mesh (defaults to 0)

        Returns:
            list(float): a list of float at each time step (according to self.solution["t"])
        """

        try:
            assert self.solve_done
        except AssertionError as err:
            logging.error("The system has not been solved yet")
            raise err
        
        try:
            assert order in [0,1,2]
        except AssertionError as err:
            logging.error(f"{err}. Unknown order {order} (0,1 or 2)")
            raise ValueError

        values = []
        for k, t in enumerate(self.solution["t"]):
            self.gf_model.set_time(t)
            # Set the solution to z^n or 0.5*(z^{n+1}+z^n)
            if CN:
                if k==0:
                    self.gf_model.to_variables(self.solution["z"][k])
                else:
                    smoothed = 0.5*(self.solution["z"][k]+self.solution["z"][k-1])
                    self.gf_model.to_variables(smoothed)
            else:
                self.gf_model.to_variables(self.solution["z"][k])
            values.append(
                gf.asm_generic(
                    self.domain._int_method[mesh_id],
                    order,
                    expression,
                    region,
                    self.gf_model,
                )
            )

        return values

    def display(self, verbose=2):
        """A method giving infos about the dphs

        Args:
            - verbose (int): the level of verbosity (defaults to 2)
        """

        # TODO: improve presentation.

        if verbose > 0:
            self.gf_model.variable_list()
            self.gf_model.brick_list()

        if verbose > 1:
            self.domain.display()

            print(self.states)
            print(self.costates)
            print(self.ports)
            print(self.bricks)
            print(self.controls)

            for term in self.hamiltonian.get_terms():
                print(term)

            for _, port in self.ports.items():
                print(port)
