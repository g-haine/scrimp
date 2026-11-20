#!/usr/bin/env python3
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Lean native GetFEM 5.5 Nedelec-II simplex smoke (experimental source patch).

Only triangles/tetrahedra and second-kind orders k=1,2. This is NOT
an upstream certification, and it must run inside the disposable patched
Docker container, not the qualified scrimp-mpi:local image directly.
"""
from __future__ import annotations

import numpy as np
import getfem as gf


def mesh_2d():
    mesh = gf.Mesh("empty", 2)
    geo = gf.GeoTrans("GT_PK(2,1)")
    # Two triangles with the shared (0,0)--(1,0) edge reversed.
    mesh.add_convex(geo, np.array([[0., 1., 0.], [0., 0., 1.]]))
    mesh.add_convex(geo, np.array([[1., 0., 0.], [0., 0., -1.]]))
    return mesh


def mesh_3d():
    mesh = gf.Mesh("empty", 3)
    geo = gf.GeoTrans("GT_PK(3,1)")
    # Same common face, opposite tetrahedron vertices and orientation.
    mesh.add_convex(geo, np.array([
        [0., 1., 0., 0.],
        [0., 0., 1., 0.],
        [0., 0., 0., 1.],
    ]))
    mesh.add_convex(geo, np.array([
        [0., 0., 1., 0.],
        [0., 1., 0., 0.],
        [0., 0., 0., -1.],
    ]))
    return mesh


def expected_local(dim, k):
    return (k + 1)*(k + 2) if dim == 2 else (k + 1)*(k + 2)*(k + 3)//2


def interpolate(mf, coefs, points):
    values = np.asarray(gf.compute(mf, coefs, "interpolate on", points),
                        dtype=float)
    if values.shape == (points.shape[1], points.shape[0]):
        values = values.T
    assert values.shape == points.shape, (
        f"Interpolation shape {values.shape}, expected {points.shape}"
    )
    assert np.all(np.isfinite(values))
    return values


def test_pair(dim: int, k: int):
    mesh = mesh_2d() if dim == 2 else mesh_3d()
    mf = gf.MeshFem(mesh, dim)
    mf.set_fem(gf.Fem(f"FEM_NEDELEC2({dim},{k})"))
    local = expected_local(dim, k)
    # Shared-edge k+1 tangential DoFs in 2D. On the shared
    # tetrahedron face, k<=2 contributes 3(k+1) edge DoFs
    # and 3(k-1) intrinsic face-tangential DoFs.
    shared = k+1 if dim == 2 else 3*(k+1) + 3*(k-1)
    expected = 2*local - shared
    n = int(mf.nbdof())
    print(f"N2 dim={dim},k={k}: global DOFs={n} expected={expected}",
          flush=True)
    assert n == expected, f"DoF connectivity failure dim={dim} k={k}"

    mim = gf.MeshIm(mesh, 2*(k+2))
    M = np.asarray(gf.asm("mass matrix", mim, mf).full(), dtype=float)
    assert M.shape == (n, n), (M.shape, n)
    symmetry = float(np.linalg.norm(M-M.T)/max(1.,np.linalg.norm(M)))
    smallest = float(np.linalg.eigvalsh((M+M.T)/2).min())
    print(f"  mass: symmetric_err={symmetry:.2e}, "
          f"lambda_min={smallest:.3e}", flush=True)
    assert symmetry < 1e-10 and smallest > 0.

    # Interpolate from inside each neighboring cell, not on the
    # geometric interface where element selection can be ambiguous.
    e = 1e-7
    x = np.array([.14, .29, .41, .57, .73])
    if dim == 2:
        pts = np.vstack([np.r_[x, x],
                         np.r_[np.full(5,e), np.full(5,-e)]])
        tangent_indices = (0,)
    else:
        y = np.array([.18, .24, .16, .21, .11])
        pts = np.vstack([np.r_[x, x], np.r_[y, y],
                         np.r_[np.full(5,e), np.full(5,-e)]])
        tangent_indices = (0,1)
    idx = np.arange(n, dtype=float)
    for j, coeffs in enumerate((
        np.sin(.61*idx)+.4*np.cos(.77*idx),
        np.cos(.47*idx)-.3*np.sin(.33*idx),
    )):
        val = interpolate(mf, coeffs, pts)
        tangential = val[list(tangent_indices)]
        jumps = tangential[:,:5]-tangential[:,5:]
        err = float(np.max(np.abs(jumps)))
        relative = err/max(1.,float(np.max(np.abs(tangential))))
        print(f"  tangential trace field={j}: abs={err:.3e}, "
              f"relative={relative:.3e}", flush=True)
        assert relative < 1e-4, (
            f"Nonconforming tangential trace dim={dim},k={k}: {relative}"
        )
    print(f"PASS: N2 dim={dim},k={k}, "
          "two-cell linkage, mass and tangential trace", flush=True)


def main():
    print("GetFEM package:", gf.__file__, flush=True)
    print("GetFEM MPI level:", gf.util_mpi_parallelism_level(), flush=True)
    assert int(gf.util_mpi_parallelism_level()) == 2
    for dim in (2,3):
        for k in (1,2):
            descriptor = f"FEM_NEDELEC2({dim},{k})"
            n = int(gf.Fem(descriptor).nbdof())
            expected = expected_local(dim,k)
            print(f"{descriptor}: local DOFs={n}, expected={expected}",
                  flush=True)
            assert n == expected
        for k in (1,2):
            test_pair(dim,k)
    print("PASS: experimental simplicial Nedelec-II low-order native smoke.")
    print("Note: numerical smoke only; not Maxwell or MPI rank-parity.")


if __name__ == "__main__":
    main()
