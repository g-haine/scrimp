# SCRIMP v1.3.0 — GetFEM 5.5 finite-element family audit

**Status:** source-level audit (19 November 2026); **native C++/Python and MPI qualification not yet performed**. Do not mistake constructor acceptance for certified H(div)/H(curl) conformity.

## Immutable source reference

- Official GetFEM Savannah repository: https://git.savannah.nongnu.org/git/getfem.git
- GitHub mirror: https://github.com/getfem/getfem
- Annotated tag: `v5.5` (8 June 2026)
- **Dereferenced release commit: `c242f7fe0989e9a90846c63eaccc6b495d774724`**
- Principal source: `src/getfem_fem.cc` in that commit; see RTk lines 1783–1953, RT0Q lines 2064–2158, BDMk lines 2161–2351, Nédélec lines 2354–2464, registry lines 4639–4644.
- Reference grids: `src/bgeot_convex_ref.cc`, `src/getfem/bgeot_convex_ref.h` and `doc/sphinx/source/project/femdesc.rst`.
- User documentation: https://getfem.org/userdoc/appendixA.html
- Mathematical cross-check: https://defelement.org/

## Family / geometry availability (GetFEM 5.5)

Here `k` is the conventional family subdegree; lowest-order first-kind Nédélec has `k=0` despite linear shape functions.

| Family | Triangle (2D) | Tetrahedron (3D) | Quadrilateral (2D) | Hexahedron (3D) | Native status |
|---|---|---|---|---|---|
| Raviart–Thomas `RT_k` | `FEM_RTK(2,k)`, k >= 0 | `FEM_RTK(3,k)`, k >= 0 | `FEM_RT0Q(2)` only k=0 | `FEM_RT0Q(3)` only k=0 | Parameterized simplex construction present, not yet runtime-qualified |
| Brezzi–Douglas–Marini `BDM_k` | `FEM_BDMK(2,k)`, k >= 1 | `FEM_BDMK(3,k)` **declared but suspicious/broken** | Not provided | Not provided | k>2 needs `GMM_USES_LAPACK`; 3D facet-cardinality discrepancy |
| Nédélec first kind `N1_k` | `FEM_NEDELEC(2)` k=0 only | `FEM_NEDELEC(3)` k=0 only | Not provided | Not provided | Whitney lowest-order simplicial construction only |
| Nédélec second kind `N2_k` | No registered constructor | No registered constructor | No registered constructor | No registered constructor | Missing |

`FEM_RT0(n)` redirects to `FEM_RTK(n,0)`. The RTK/BDMK constructors accept integers `2 <= n < 100`, `0 <= k < 100` (RT), and `1 <= k < 100` (BDM). This is **parameter validation**, not an arbitrary-order numerical guarantee. The registration contains `NEDELEC` but no higher-order or second-kind Nédélec suffix.

## Reference dimensions — one cell, not assembled mesh DoFs

| Space | Triangle | Tetrahedron | Order |
|---|---|---|---|
| `RT_k = [P_k]^d + x \\widetilde{P}_k` | `(k+1)(k+3)` | `(k+1)(k+2)(k+4)/2` | k >= 0 |
| `BDM_k = [P_k]^d` | `(k+1)(k+2)` | `(k+1)(k+2)(k+3)/2` | k >= 1 |
| `N1_k = [P_k]^d \\oplus Z_{k+1}, Z_{k+1}={p homogeneous of degree k+1: x.p=0}` | `(k+1)(k+3)` | `(k+1)(k+3)(k+4)/2` | k >= 0 |
| `N2_k = [P_k]^d` | `(k+1)(k+2)` | `(k+1)(k+2)(k+3)/2` | k >= 1 |

Examples: RT0 2D/3D = 3/4; N1_0 2D/3D = 3/6; BDM1 2D/3D = 6/12; N2_1 2D/3D = 6/12. For Maxwell 3D, distinguish polynomial subdegree from the degree of the shape functions.

## Source findings requiring native checks

1. **RT_k on simplices exists for any descriptor order accepted by the constructor.** It uses a polynomial basis and nodal normal/facet DOFs on a `k+n` simplex Lagrange lattice, plus interior component evaluations. Its DOF functionals are not the canonical integral moments; confirm polynomial trace completeness, physical facet matching and exact Piola behavior, particularly at high k.
2. **BDM_k in 3D has a likely defect.** The constructor uses `simplex_of_reference(n, k+2)` and creates facet-normal DOFs only at lattice points belonging to exactly one facet. For a tetrahedron with k=1, that produces just 4 facet points, while BDM1 has 12 DOFs and adds no interior rows. At k=2 it creates 12 facet rows plus 4 divergence rows = 16, while BDM2 has 30 DoFs. The `ipoint == nddl` assertion therefore appears unreachable in 3D k=1,2. **Static prediction; run isolated native constructor tests before classifying a runtime failure.** At k>=3, an SVD may complete the square matrix without restoring the missing normal facet moments; check H(div) conformity independently.
3. **BDM_k, k>2 uses SVD conditional on `GMM_USES_LAPACK`.** The code uses absolute tolerances 1e-14 and 1e-16, a potential conditioning/rank weakness at high order. Do not characterize k=3..5 as validated.
4. **N1_0 uses standard barycentric Whitney 1-form shapes** `lambda_i grad(lambda_j)-lambda_j grad(lambda_i)`, one edge DOF per edge and a geometry-dependent `mat_trans` matrix. Higher-order edge/face/interior moment classes and the second-kind Nédélec family are absent from this registry.
5. **Orientation/transformation audit:** RT/BDM `mat_trans` uses a fixed direction `(pi,pi^2,pi^3,...)` to choose signs of transformed normals, and Nédélec similarly orients tangents. Verify mirrored/reordered cells, shared facets/edges, and both Piola pullbacks. Source loops in RT/BDM use `j <= nb_dof(0)`: inspect `faces_of_dof` bounds before concluding whether this is an off-by-one defect.
6. **Geometry restriction:** `FEM_RT0Q` is an independently implemented low-order tensor-product H(div) element. It does not imply higher-order RT on quads/hexes. Do not attempt simplex descriptors on tensor cells.

For comparison, conventional commuting pullbacks are:
- H(div), contravariant: `v(x) = J vhat(xhat) / det(J)` for orientation-preserving maps (signed determinant convention needs explicit orientation tests).
- H(curl), covariant: `e(x) = J^{-T} ehat(xhat)`.

## Native audit gates before any implementation

Use **GetFEM v5.5 original sources** and the local SCRIMP Docker/GetFEM build; capture compiler flags, `GMM_USES_LAPACK`, Python/C++ versions, and the exact GetFEM commit. Run constructor probes in separate processes (an upstream assertion can abort the process).

- **A1, descriptors and per-cell dimensions:** in 2D and 3D test RT k=0..5, BDM k=1..5, N1 k=0; attempt N1 higher and N2 only as expected-negative descriptors. Test `FEM_RT0Q(2/3)` separately. Compare `getfem.Fem(...).nbdof()` with reference formulas, record constructor errors rather than skipping them.
- **A2, local interpolation/unisolvance:** reference shape polynomial degree, Kronecker/dual functionals, matrix rank, shape evaluation at random interior points; track conditioning with increasing k.
- **A3, conformity:** two nonidentically oriented triangles and tetrahedra; verify continuous shared normal traces (Hdiv) or tangential traces (Hcurl), including explicit change of vertex numbering and physical Jacobians.
- **A4, commuting/compatibility:** polynomial reproduction, `div(RT_k) = P_k` elementwise, `div(BDM_k) = P_{k-1}`; for compatible Hcurl spaces check grad inclusions and discrete curl/grad, div/curl identities. Separate polynomial-space statements from global mesh test results.
- **A5, numerical convergence:** manufactured 2D/3D Hdiv mixed and Hcurl problems, h refinement and at least two polynomial degrees, with expected error slopes and trace/port errors.
- **A6, MPI parity:** identical global meshes with 1/2/4 ranks; compare global dimensions, assembled mass and differential operators, sparse coefficients, full vectors, discrete energy/power residuals. Preserve SCRIMP point-0 regression baselines.
- **A7, geometry and build permutations:** affine/reflected simplexes, curved maps where supported, quad/hexa RT0Q; compare GetFEM with/without LAPACK, and note any unsupported combinations.

**Priority:** reproduce the predicted BDM1/2 tetrahedral failure first, then qualify existing RT_k before writing any replacement code. Do not patch original GetFEM in the production image before an isolated baseline is saved.

## Repository strategy and upstream contribution route

- SCRIMP point 0 is frozen at `c09a234048ec60cea928f83dfc5269b381d3581e`; `release/1.3.0-preparation` is branched directly from it. Do **not** merge SCRIMP to `master` until final release gates are satisfied.
- **Keep GetFEM source changes in a different repository.** `getfem/getfem` on GitHub is a *mirror*. Official instructions say Savannah is the main repository and recommend a `devel-<login>-<subject>` branch, project membership/SSH key, and maintainer merge request (https://getfem.org/project/contribute.html).
- Until contributor access is established, work on a separate GetFEM clone or GitHub fork pinned to `v5.5`, with individually reviewable patches: (i) existing BDM correctness / orientation fixes, (ii) Nédélec first-kind k>=1, (iii) Nédélec second kind k>=1, plus regression tests.
- Pin the verified **commit + patch SHA** in Docker only after native tests; record whether the patch is locally vendored or upstream accepted.
- Avoid coupling the GetFEM implementation to SCRIMP internals; add small SCRIMP smoke/integration tests only.

## Sources

- https://github.com/getfem/getfem/blob/v5.5/src/getfem_fem.cc
- https://github.com/getfem/getfem/blob/v5.5/src/bgeot_convex_ref.cc
- https://github.com/getfem/getfem/blob/v5.5/doc/sphinx/source/project/femdesc.rst
- https://getfem.org/userdoc/appendixA.html
- https://defelement.org/elements/raviart-thomas.html
- https://defelement.org/elements/brezzi-douglas-marini.html
- https://defelement.org/elements/nedelec1.html
- https://defelement.org/elements/nedelec2.html
- https://getfem.org/project/contribute.html

## Native constructor probe — container qualification (19 November 2026)

The historical host-conda execution on 2026-11-19T10:38Z loaded the **old** GetFEM from `/opt/anaconda3/envs/scrimp`. Its 5/14 result does **not** qualify GetFEM 5.5. Keep the JSON as a labeled historical non-qualifying probe.

The subsequent **real Docker** execution loaded Python and GetFEM from `/opt/conda/envs/scrimp`, reported `docker_marker_present=true`, and `getfem.util_mpi_parallelism_level() == 2`. It constructed **12 of 14** descriptors with matching local dimensions:

| Family / dimension | k range tested | Constructed / attempted | Conclusion |
|---|---:|---:|---|
| RT simplex 2D | 0–2 | 3/3 | Local construction and dimension PASS |
| RT simplex 3D | 0–2 | 3/3 | Local construction and dimension PASS |
| BDM simplex 2D | 1–2 | 2/2 | Local construction and dimension PASS |
| **BDM simplex 3D** | **1–2** | **0/2** | Constructor ERROR; detailed exception not yet captured |
| RT0Q 2D and 3D | 0 | 2/2 | Local construction and dimension PASS |
| Nédélec first kind 2D and 3D | 0 | 2/2 | Local construction and dimension PASS |

The first Docker script wrote the JSON **after** all probes to `sandbox/outputs/getfem55_constructor_audit_docker.json` on the bind-mounted host tree. This failed with `PermissionError: [Errno 13]` because the non-root image user lacked write permission to the destination; the child error messages for BDM were therefore not preserved in JSON. **This is a report-output permission issue, not a GetFEM constructor failure.**

### Corrected reporting without changing host permissions

The audit script now supports `--output -`: diagnostics go to **stderr**, while complete machine-readable JSON goes to **stdout**. Shell redirection on the **host** writes it to a host-owned path, so no Docker user needs write access to `sandbox/outputs`. Unwritable file output also falls back to printing the JSON between markers instead of raising a traceback. All failed descriptors display their full exception text.

From the updated `release/1.3.0-preparation` worktree, test **only the failed constructors**:

```bash
cd ~/Simulations/scrimp
git pull --ff-only
sudo docker run --rm \
  -v "$PWD:/workspace/scrimp:ro" \
  -w /workspace/scrimp \
  scrimp-mpi:local \
  python -u sandbox/audit_getfem55_element_catalog.py \
  --max-order 2 \
  --only 'FEM_BDMK(3,1)' \
  --only 'FEM_BDMK(3,2)' \
  --output - \
  > "$HOME/getfem55_bdm3_docker.json"
```

Once the BDM error traces are inspected, broaden to `--max-order 5` as resources permit; capture to a **different** JSON filename. Record the Docker image ID (`sudo docker image inspect scrimp-mpi:local --format '{{.Id}}'`) and build provenance (actual GetFEM source archive checksum/version). The runtime gate proves Python/extension origin and MPI level 2, but it **cannot independently prove** the source commit: `source_commit_verified=false`.

**Scope:** local element construction only, not H(div)/H(curl) conformity, shared facets/edge orientation, Piola transformations, numerical convergence, or MPI multi-rank parity. Nothing in the point-0 solver or benchmarks was changed.
