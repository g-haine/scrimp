# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Universal finite-element L2 initial-value projection (default).

GetFEM assembles the physical mass matrix and source on the *actual*
scalar/vector/derivative finite-element basis, including filtered regions.
PETSc owns the sparse linear solve. No nodal coefficient assumptions.

Historical interpolation is explicitly selectable since non-polynomial
nodal interpolation differs from L2 projection. Distributed MPI1/2 reference IJV is qualified; the new default
COO mass transfer must be qualified natively. DAE/PDE physics is separate.
"""

import ast
import math
import os

import getfem as gf
import numpy as np
from petsc4py import PETSc

from scrimp.utils.linalg import convert_gmm_to_petsc


_SCALAR_FUNCTIONS = frozenset((
    "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
    "exp", "log", "log10", "sqrt", "abs", "sqr",
))


def python_expression_to_gwfl(expression, qdim):
    """Translate a restricted subset of existing SCRIMP Python expressions.

    Supported: real scalars, x/y/z, arithmetic, power, common
    NumPy/math scalar functions including np.sign, and literal
    Python list/tuple/np.array vectors.
    This is a strict AST translator, NOT Python eval and NOT a general
    NumPy interpreter. Unsupported constructs fail before model mutation.
    """
    if not isinstance(expression, str) or not expression.strip():
        raise ValueError("Initial FE expression must be a nonempty string")
    if not isinstance(qdim, int) or qdim < 1:
        raise ValueError("Invalid initial FE field dimension")

    def call_name(f):
        if isinstance(f, ast.Name):
            return f.id
        if (isinstance(f, ast.Attribute) and
                isinstance(f.value, ast.Name) and
                f.value.id in ("np", "numpy", "math")):
            return f.attr
        raise ValueError("Unsupported initial FE function")

    def value(n):
        if isinstance(n, ast.Constant):
            if isinstance(n.value, bool) or not isinstance(
                n.value, (int, float)
            ) or not math.isfinite(float(n.value)):
                raise ValueError("Initial FE constants must be finite real numbers")
            return repr(n.value)
        if isinstance(n, ast.Name):
            if n.id in ("x", "y", "z"):
                return f"X({('x', 'y', 'z').index(n.id) + 1})"
            if n.id == "pi":
                return repr(math.pi)
            if n.id == "e":
                return repr(math.e)
            raise ValueError(f"Unsupported initial FE identifier: {n.id}")
        if isinstance(n, ast.Attribute):
            if (isinstance(n.value, ast.Name) and
                    n.value.id in ("np", "numpy", "math") and
                    n.attr in ("pi", "e")):
                return repr(getattr(math, n.attr))
            raise ValueError("Unsupported initial FE attribute")
        if isinstance(n, ast.UnaryOp):
            if isinstance(n.op, (ast.UAdd, ast.USub)):
                op = "+" if isinstance(n.op, ast.UAdd) else "-"
                return f"({op}({value(n.operand)}))"
            raise ValueError("Unsupported initial FE unary operator")
        if isinstance(n, ast.BinOp):
            if isinstance(n.op, ast.Pow):
                return f"pow({value(n.left)},{value(n.right)})"
            operation = {
                ast.Add: "+", ast.Sub: "-", ast.Mult: "*",
                ast.Div: "/",
            }.get(type(n.op))
            if operation is None:
                raise ValueError("Unsupported initial FE binary operator")
            return f"(({value(n.left)}){operation}({value(n.right)}))"
        if isinstance(n, (ast.List, ast.Tuple)):
            return "[" + ";".join(value(el) for el in n.elts) + "]"
        if isinstance(n, ast.Call):
            name = call_name(n.func)
            if n.keywords:
                raise ValueError("Initial FE function keyword arguments unsupported")
            if name == "array":
                if len(n.args) != 1 or not isinstance(n.args[0], (
                    ast.List, ast.Tuple
                )):
                    raise ValueError("np.array requires a literal vector")
                return value(n.args[0])
            if name == "sign" and len(n.args) == 1:
                # GetFEM's native sign(0) is +1, whereas np.sign(0)
                # is 0. Preserve the published shallow-water initial
                # condition *exactly*, including at the jump x=0.5.
                # GWFL Heaviside(0)=1, so H(t)-H(-t) is -1/0/+1.
                argument = value(n.args[0])
                return (
                    f"(Heaviside({argument})"
                    f"-Heaviside(-({argument})))"
                )
            if name == "power":
                name = "pow"
            if name == "pow" and len(n.args) == 2:
                return f"pow({value(n.args[0])},{value(n.args[1])})"
            if name in _SCALAR_FUNCTIONS and (
                len(n.args) == (2 if name == "atan2" else 1)
            ):
                return f"{name}(" + ",".join(value(x) for x in n.args) + ")"
            raise ValueError(f"Unsupported initial FE function: {name}")
        raise ValueError(
            f"Unsupported Python syntax in initial FE expression: "
            f"{type(n).__name__}"
        )

    tree = ast.parse(expression, mode="eval").body
    if isinstance(tree, (ast.List, ast.Tuple)) or (
        isinstance(tree, ast.Call) and call_name(tree.func) == "array"
    ):
        if isinstance(tree, (ast.List, ast.Tuple)):
            elements = tree.elts
        else:
            if (len(tree.args) != 1 or tree.keywords or
                    not isinstance(tree.args[0], (ast.List, ast.Tuple))):
                raise ValueError(
                    "np.array initial FE field must contain a literal vector"
                )
            elements = tree.args[0].elts
        if len(elements) != qdim:
            raise ValueError(
                f"Initial vector contains {len(elements)} components, "
                f"expected {qdim}"
            )
        return value(elements[0]) if qdim == 1 else value(tree)
    if qdim != 1:
        raise ValueError(
            f"Expected an explicit {qdim}-component vector initial value"
        )
    return value(tree)


def real_getfem_projection_coefficients(values, *, tolerance=1e-12):
    """Bridge PETSc scalar coefficients to a real GetFEM model safely.

    PETSc may use complex ScalarType even when the GetFEM model and all
    physical FE data are strictly real. GetFEM's `Model.set_variable`
    accepts a DOUBLE REAL array, not complex128 (even with zero
    imaginary part). Reject significant imaginary components instead
    of discarding them silently, and always return contiguous float64.
    """
    a = np.asarray(values)
    if a.ndim != 1 or not a.size or not np.issubdtype(a.dtype, np.number):
        raise ValueError("GetFEM initial coefficients must be a nonempty 1D numeric array")
    if not np.all(np.isfinite(a)):
        raise ValueError("Non-finite PETSc initial projection coefficient")
    if (not np.isfinite(tolerance) or tolerance <= 0):
        raise ValueError("Invalid PETSc-to-GetFEM real conversion tolerance")
    if np.iscomplexobj(a):
        re = np.asarray(a.real)
        im = np.asarray(a.imag)
        bound = tolerance * max(1.0, float(np.max(np.abs(re))))
        max_im = float(np.max(np.abs(im)))
        if max_im > bound:
            raise ValueError(
                "PETSc FE projection has nonnegligible imaginary data "
                f"({max_im:g} > {bound:g}); real GetFEM cannot store it"
            )
        a = re
    return np.ascontiguousarray(a, dtype=np.float64)


def scale_getfem_projection_rhs(values, rhs_scale):
    """Apply SCRIMP's *multiplicative* calibrated GetFEM RHS factor.

    _getfem_mpi_rhs_scale() returns 1/P when an MPI-enabled GetFEM
    assembly already duplicated a fully-reduced RHS P times. A
    division here would amplify the duplication rather than fix it.
    """
    scale = float(rhs_scale)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Invalid GetFEM projection RHS normalization")
    rhs = np.asarray(values, dtype=np.float64)
    if rhs.ndim != 1 or not np.all(np.isfinite(rhs)):
        raise ValueError("GetFEM projection RHS must be a finite 1D vector")
    result = np.ascontiguousarray(rhs * scale)
    if not np.all(np.isfinite(result)):
        raise ValueError("Non-finite normalized GetFEM projection RHS")
    return result


def fe_projection_transfer_policy(mpi_size, matrix_mode):
    """Select the same linear-mass GetFEM/PETSc bridge as DPHS.

    With distributed GetFEM assembly, linear mass matrices use the
    standard COO-eligible path of convert_gmm_to_petsc. PETSc assembles
    the resulting AIJ (CSR-like) matrix. The existing bridge retains
    its explicit legacy IJV override for diagnostics; this is not a
    separate numerical FE initialization method.

    MPI>1 cannot treat rank-local GetFEM contributions as a complete
    replicated global matrix. In serial, allow the historical
    replicated mode for compatibility with existing applications.
    """
    if isinstance(mpi_size, bool) or int(mpi_size) < 1:
        raise ValueError("Invalid PETSc MPI size for FE initialization")
    mode = str(matrix_mode).strip().lower()
    if mode not in ("distributed", "replicated"):
        raise ValueError("Unknown GetFEM matrix mode for FE initialization")
    if int(mpi_size) > 1 and mode != "distributed":
        raise ValueError(
            "L2 FE initialization with MPI > 1 requires "
            "SCRIMP_GETFEM_MATRIX_MODE=distributed"
        )
    return mode == "distributed"

def project_fe_initial_value(
    name, expression, fem, region, mim, *,
    rhs_scale=1.0, relative_tolerance=1e-9,
):
    """Compute L2 projection into the *exact* GetFEM FE coefficient basis.

    All ranks participate. With MPI>1 the matrix bridge requires the
    established distributed GetFEM assembly semantics. Linear mass
    matrices use the same COO-eligible PETSc bridge as DPHS. One global RHS is inserted on PETSc
    owned rows and the coefficient array is gathered for the original
    GetFEM model interface.
    """
    coo_eligible = fe_projection_transfer_policy(
        PETSc.COMM_WORLD.getSize(),
        os.environ.get("SCRIMP_GETFEM_MATRIX_MODE", "distributed"),
    )
    qdim = int(fem.qdim())
    source = python_expression_to_gwfl(expression, qdim)
    if not np.isfinite(float(rhs_scale)) or float(rhs_scale) <= 0:
        raise ValueError("Invalid GetFEM MPI source calibration")
    if (not np.isfinite(float(relative_tolerance)) or
            not 0 < float(relative_tolerance) < 1):
        raise ValueError("Invalid FE projection residual tolerance")

    model = gf.Model("real")
    if region is None:
        model.add_fem_variable("initial_projection", fem)
        region_id = -1
    else:
        region_id = int(region)
        model.add_filtered_fem_variable("initial_projection", fem, region_id)

    scalar = qdim == 1
    mass_form = (
        "initial_projection*Test_initial_projection" if scalar
        else "initial_projection.Test_initial_projection"
    )
    source_form = (
        f"({source})*Test_initial_projection" if scalar
        else f"({source}).Test_initial_projection"
    )
    model.add_linear_term(mim, mass_form, region_id)
    model.add_source_term(mim, source_form, region_id)
    model.assembly()
    n = int(model.nbdof())
    if n <= 0:
        raise ValueError(f"No FE coefficients for initial variable {name}")
    if region is not None:
        region_dofs = np.asarray(fem.basic_dof_on_region(region_id)).size
        if region_dofs != n:
            raise ValueError(
                f"Filtered FE DOF count mismatch for {name}: "
                f"GetFEM model {n} vs mesh fem {region_dofs}"
            )
    gf_rhs = scale_getfem_projection_rhs(
        np.asarray(model.rhs(), dtype=float).reshape(-1), rhs_scale
    )
    if gf_rhs.size != n or not np.all(np.isfinite(gf_rhs)):
        raise ValueError(f"Invalid GetFEM projection right-hand side for {name}")

    A = b = x = ksp = residual = Ax = scatter = gathered = None
    try:
        A = PETSc.Mat().createAIJ(size=(n, n), comm=PETSc.COMM_WORLD)
        convert_gmm_to_petsc(
            model.tangent_matrix(), A, comm=PETSc.COMM_WORLD,
            distributed_coo_eligible=coo_eligible,
        )
        b = A.createVecLeft()
        x = A.createVecRight()
        lo, hi = b.getOwnershipRange()
        if hi > lo:
            b.setValues(
                np.arange(lo, hi, dtype=PETSc.IntType),
                gf_rhs[lo:hi],
                addv=PETSc.InsertMode.INSERT_VALUES,
            )
        b.assemble()
        ksp = PETSc.KSP().create(comm=PETSc.COMM_WORLD)
        ksp.setOptionsPrefix("scrimp_fe_initial_")
        ksp.setOperators(A)
        ksp.setType(PETSc.KSP.Type.PREONLY)
        pc = ksp.getPC()
        pc.setType(PETSc.PC.Type.LU)
        pc.setFactorSolverType("mumps")
        ksp.setFromOptions()
        ksp.solve(b, x)
        if int(ksp.getConvergedReason()) <= 0:
            raise RuntimeError(
                f"PETSc FE projection KSP failed: {ksp.getConvergedReason()}"
            )
        Ax = b.duplicate()
        A.mult(x, Ax)
        residual = b.copy()
        residual.axpy(-1.0, Ax)
        rel = float(residual.norm()) / max(float(b.norm()), 1e-30)
        if not np.isfinite(rel) or rel > relative_tolerance:
            raise RuntimeError(
                f"FE initial projection residual {rel:.6g} exceeds "
                f"{relative_tolerance:.6g} for {name}"
            )
        scatter, gathered = PETSc.Scatter().toAll(x)
        scatter.scatter(x, gathered)
        coefficients = real_getfem_projection_coefficients(
            gathered.getArray()
        )
        if coefficients.size != n:
            raise RuntimeError("Invalid gathered FE projection coefficients")
        return coefficients
    finally:
        for object_ in (gathered, scatter, residual, Ax, ksp, x, b, A):
            if object_ is not None:
                object_.destroy()


def filtered_dirichlet_dofs(mesh_fem, volume_region, boundary_region,
                            *, coefficient_count):
    """Map a conforming scalar FE boundary to its filtered model indices.

    A filtered GetFEM variable has a different coefficient ordering from
    the global MeshFem. Use the actual GetFEM basic_dof_on_region order,
    not mesh points, component interleaving or a guessed CG node index.
    The calling model must guarantee the boundary FEM trace is essential
    (e.g. a scalar CG1). Do not use this for H(div)/H(curl) normal data.
    """
    domain = np.asarray(
        mesh_fem.basic_dof_on_region(volume_region), dtype=np.int64
    ).reshape(-1)
    boundary = np.asarray(
        mesh_fem.basic_dof_on_region(boundary_region), dtype=np.int64
    ).reshape(-1)
    if (domain.size != coefficient_count or not boundary.size
            or np.unique(domain).size != domain.size
            or np.unique(boundary).size != boundary.size
            or np.setdiff1d(boundary, domain).size):
        raise ValueError("Filtered GetFEM domain/boundary DOF partition invalid")
    fixed = np.flatnonzero(np.isin(domain, boundary))
    free = np.flatnonzero(~np.isin(domain, boundary))
    if not fixed.size or not free.size:
        raise ValueError("Essential Dirichlet boundary must have fixed/free DOFs")
    return fixed, free


def constrained_l2_zero_boundary_correction(
    system, name, initial, fixed, free, *, volume_region,
    factor_solver="mumps",
):
    """Mass-orthogonal correction of an already L2-projected CG field.

    For original coefficients c0 and exact physical FE mass M, solve
    M_FF (c-c0)_F = M_FB (c0)_B, subject to c_B=0.
    This minimizes ||c-c0||_M with an exact essential zero boundary.
    GetFEM's mass, sparse COO PETSc bridge and MUMPS are used, with a
    native stationarity check. This is the independently qualified
    Heat-Wave constrained-L2 initial-data correction, not pseudo-DAE.
    """
    if PETSc.COMM_WORLD.getSize() != 1:
        raise NotImplementedError(
            "Constrained-L2 boundary correction requires MPI1; "
            "use the already MPI1/2-qualified clamp for MPI>1"
        )
    fem = system.ports[name].get_fem()
    mim = system.domain._int_method[system.ports[name].get_mesh_id()]
    model = gf.Model("real")
    model.add_filtered_fem_variable("boundary_corrected", fem, volume_region)
    model.add_linear_term(
        mim, "boundary_corrected*Test_boundary_corrected", volume_region
    )
    model.assembly()
    n = int(initial.size)
    if int(model.nbdof()) != n:
        raise ValueError("Wrong physical filtered GetFEM mass dimension")

    mat = PETSc.Mat().createAIJ(size=(n, n), comm=PETSc.COMM_WORLD)
    x_boundary = m_boundary = subset = m_free = rhs = delta = ksp = None
    correction_vec = mass_correction = None
    try:
        convert_gmm_to_petsc(
            model.tangent_matrix(), mat, comm=PETSc.COMM_WORLD,
            distributed_coo_eligible=True,
        )
        x_boundary = mat.createVecRight()
        x_boundary.set(0.)
        x_boundary.setValues(
            fixed.astype(PETSc.IntType), initial[fixed],
            addv=PETSc.InsertMode.INSERT_VALUES,
        )
        x_boundary.assemble()
        m_boundary = mat.createVecLeft()
        mat.mult(x_boundary, m_boundary)
        subset = PETSc.IS().createGeneral(
            free.astype(PETSc.IntType), comm=PETSc.COMM_WORLD
        )
        m_free = mat.createSubMatrix(subset, subset)
        rhs = m_free.createVecLeft()
        rhs.setValues(
            np.arange(len(free), dtype=PETSc.IntType),
            m_boundary.getValues(free.astype(PETSc.IntType)),
            addv=PETSc.InsertMode.INSERT_VALUES,
        )
        rhs.assemble()
        delta = m_free.createVecRight()
        ksp = PETSc.KSP().create(comm=PETSc.COMM_WORLD)
        ksp.setOptionsPrefix("scrimp_initial_dirichlet_l2_")
        ksp.setOperators(m_free)
        ksp.setType(PETSc.KSP.Type.PREONLY)
        ksp.getPC().setType(PETSc.PC.Type.LU)
        ksp.getPC().setFactorSolverType(factor_solver)
        ksp.setFromOptions()
        ksp.solve(rhs, delta)
        if int(ksp.getConvergedReason()) <= 0:
            raise RuntimeError("Constrained-L2 boundary mass KSP did not converge")
        delta_free = real_getfem_projection_coefficients(
            np.asarray(delta.getArray(readonly=True))
        )
        corrected = initial.copy()
        corrected[fixed] = 0.
        corrected[free] += delta_free
        if not np.all(np.isfinite(corrected)) or np.any(corrected[fixed]):
            raise RuntimeError("Boundary-corrected FE field violates zero trace")
        correction_vec = mat.createVecRight()
        correction_vec.setValues(
            np.arange(n, dtype=PETSc.IntType), corrected-initial,
            addv=PETSc.InsertMode.INSERT_VALUES,
        )
        correction_vec.assemble()
        mass_correction = mat.createVecLeft()
        mat.mult(correction_vec, mass_correction)
        derivative = real_getfem_projection_coefficients(
            np.asarray(mass_correction.getArray(readonly=True))
        )
        stationarity = float(np.linalg.norm(derivative[free]) /
                             max(float(np.linalg.norm(derivative)), 1e-30))
        mass_squared = float(np.vdot(
            correction_vec.getArray(readonly=True),
            mass_correction.getArray(readonly=True),
        ).real)
        if (stationarity > 1e-8 or not np.isfinite(mass_squared)
                or mass_squared < -1e-12):
            raise RuntimeError(
                "Constrained-L2 correction fails physical mass stationarity"
            )
        return corrected, {
            "mass_norm_correction": float(np.sqrt(max(mass_squared, 0.))),
            "relative_stationarity_on_free_dofs": stationarity,
            "interior_mass_solve_KSP_reason": int(ksp.getConvergedReason()),
        }
    finally:
        for item in (
            mass_correction, correction_vec, ksp, delta, rhs,
            m_free, subset, m_boundary, x_boundary, mat,
        ):
            if item is not None:
                item.destroy()


def enforce_zero_dirichlet_initial_values(
    system, name, *, volume_region, boundary_region,
    method="clamp", factor_solver="mumps",
):
    """Correct only the already-projected CG wave initial coefficients.

    First project the prescribed analytic field via the usual
    DPHS.set_initial_value(..., method="l2"). Then ensure its *exact
    discrete* trace vanishes. 'clamp' preserves interior coefficients;
    'constrained_l2' changes interior coefficients minimally in the
    physical GetFEM mass norm. Neither invokes the DAE initializer or
    alters the heat/stress state. This belongs to the *model's*
    physical initial data, NOT an automatic relaxation of the DAE.
    """
    if method not in ("clamp", "constrained_l2"):
        raise ValueError("Dirichlet initial data: use clamp or constrained_l2")
    if name not in system.ports:
        raise ValueError("Missing finite-element port for initial boundary data")
    fem = system.ports[name].get_fem()
    if int(fem.qdim()) != 1:
        raise ValueError("Only scalar essential Dirichlet traces supported")
    initial = real_getfem_projection_coefficients(
        np.asarray(system.gf_model.variable(name)).reshape(-1)
    )
    fixed, free = filtered_dirichlet_dofs(
        fem, volume_region, boundary_region, coefficient_count=len(initial),
    )
    details = {}
    if method == "clamp":
        corrected = initial.copy()
        corrected[fixed] = 0.
        if not np.array_equal(corrected[free], initial[free]):
            raise RuntimeError("Clamping changed an interior FE coefficient")
    else:
        corrected, details = constrained_l2_zero_boundary_correction(
            system, name, initial, fixed, free,
            volume_region=volume_region, factor_solver=factor_solver,
        )
    if (not np.all(np.isfinite(corrected))
            or np.any(corrected[fixed])):
        raise RuntimeError("Dirichlet initial field has nonzero boundary DOFs")
    system.set_from_vector(name, corrected)
    if not np.array_equal(
        np.asarray(system.gf_model.variable(name)).reshape(-1), corrected
    ):
        raise RuntimeError("GetFEM rejected physical compatible initial state")
    return {
        "method": method,
        "field": name,
        "volume_region": int(volume_region),
        "boundary_region": int(boundary_region),
        "filtered_DOFS": int(len(initial)),
        "boundary_DOFS": int(len(fixed)),
        "free_DOFS": int(len(free)),
        "before_boundary_l2_coefficients": float(np.linalg.norm(initial[fixed])),
        "after_boundary_l2_coefficients": float(np.linalg.norm(corrected[fixed])),
        "relative_change": float(
            np.linalg.norm(corrected-initial) /
            max(float(np.linalg.norm(initial)), 1e-30)
        ),
        "T_and_q_unmodified_by_correction": True,
        "dae_tolerances_modified": False,
        "correction": details,
    }
