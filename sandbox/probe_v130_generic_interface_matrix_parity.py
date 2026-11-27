# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Full sparse K parity: legacy Concentric helper versus generic topology.

This sandbox-only wrapper reuses the same published Heat-Wave example,
the same original weak forms, the same manually assigned interface port
and brick regions, and the previously qualified assembly probe. ONLY
the generator for side-30/side-31 face regions changes. A full sparse
matrix is exported for numerical comparison, not merely Frobenius norms.

The helper is NOT an automatic production DPHS interface compiler.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix, save_npz

from scrimp.interface import resolve_conforming_interface
import probe_v130_hw_dual_sided_assembly as original


def generic_faces(mesh):
    if set(np.asarray(mesh.regions(), dtype=np.int64).reshape(-1)) & {30, 31}:
        raise RuntimeError("Reference region labels 30/31 are occupied")
    logical = np.asarray(mesh.region(10), dtype=np.int64).copy()
    sides = resolve_conforming_interface(mesh, 2, 10, (1, 2))
    heat = sides.for_subdomain(1)
    wave = sides.for_subdomain(2)
    if heat.shape != (2, 78) or wave.shape != (2, 78):
        raise RuntimeError("Unexpected original Heat-Wave facet count")
    if set(map(tuple, heat.T)) & set(map(tuple, wave.T)):
        raise RuntimeError("Resolved oriented facets overlap")
    mesh.set_region(30, heat)
    mesh.set_region(31, wave)
    if not np.array_equal(logical, mesh.region(10)):
        raise RuntimeError("Modified original logical interface region")
    return {
        "method": "scrimp.interface.resolve_conforming_interface",
        "logical_interface_region": 10,
        "side_regions": [30, 31],
        "side_counts": [int(heat.shape[1]), int(wave.shape[1])],
        "region10_unchanged": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", choices=("reference", "generic"),
                        required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists() or args.out.with_suffix(".npz").exists():
        raise RuntimeError("Refuse to overwrite matrix evidence")
    saved_argv = sys.argv
    saved_geometry = original.new_one_sided_regions
    saved_report = original.compare_block_norms
    called = []

    def record_full_matrix(system, variant, topology, modifications, out):
        if called:
            raise RuntimeError("More than one native assembly captured")
        saved_report(system, variant, topology, modifications, out)
        n = int(system.gf_model.nbdof())
        indptr, indices, values = system.stiffness.getValuesCSR()
        K = csr_matrix(
            (np.asarray(values).copy(),
             np.asarray(indices).copy(),
             np.asarray(indptr).copy()),
            shape=(n, n),
        )
        K.sum_duplicates()
        K.sort_indices()
        save_npz(args.out.with_suffix(".npz"), K, compressed=True)
        called.append({"shape": list(K.shape), "nnz": int(K.nnz)})
        print("FULL_MATRIX_EXPORTED", json.dumps(called[-1]), flush=True)

    try:
        if args.geometry == "generic":
            original.new_one_sided_regions = generic_faces
        original.compare_block_norms = record_full_matrix
        sys.argv = [
            saved_argv[0], "--variant", "dual_sided",
            "--out", str(args.out),
        ]
        original.main()
    finally:
        sys.argv = saved_argv
        original.new_one_sided_regions = saved_geometry
        original.compare_block_norms = saved_report
    if len(called) != 1 or not args.out.with_suffix(".npz").is_file():
        raise RuntimeError("Full native matrix export did not complete")


if __name__ == "__main__":
    main()
