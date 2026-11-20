# SCRIMP v1.3.0 — Ubuntu / Docker / VS Code developer-workflow acceptance

**Date:** 2026-11-20  
**Status:** **PASS — developer workflow on the tested Ubuntu machine**  
**Branch:** `release/1.3.0-preparation`  
**Tested candidate image:** `scrimp-mpi:v1.3.0-rc1`  
**Image source commit:** `96961d647c2caf41dd77d2a7b5926e86e86c4b69`  
**Image ID / manifest-list digest reported by Docker:** `sha256:ab08d9beaa60a1685b7e8df1e69f5a085f5571ca46c2b7d8554f53ec42ace3e8`

This acceptance is based on the developer's **actual Ubuntu VS Code and Docker execution**, not a claim of automated GUI tests or a general platform release.

## Validated user path

| Acceptance check | Evidence | Result |
| --- | --- | --- |
| GitHub-driven isolated build context | Build transferred about 1.60 MB (rather than untracked simulation data) | **PASS** |
| Candidate Docker image builds | GetFEM 5.5, complex CPU PETSc/SLEPc 3.26.1, NumPy 2 and Python bindings | **PASS** |
| Native FEM and solver smoke | GetFEM MPI level 2, tested RT/BDM/Nédélec-I/II constructors, distributed MUMPS at 1 and 2 ranks | **PASS** |
| Independent project creation | `scripts/create-scrimp-project.py` creates only the documented starter files | **PASS** |
| VS Code Dev Container startup | Image `scrimp-mpi:v1.3.0-rc1`, user `mambauser`, host project bind-mounted to `/work` | **PASS** |
| Python / SCRIMP runtime | `python simulations/first_simulation.py`: 10 accepted steps, zero nonlinear failures | **PASS** |
| Persistent outputs | 11 time points with `time,q,p,H` saved to host `outputs/oscillator.csv` | **PASS** |
| Two-rank MPI startup | `mpiexec -n 2 python simulations/mpi_check.py` reports rank 1/2, rank 2/2, PASS | **PASS** |
| Python debugger | F5 breakpoint, F10 step over, Locals and Call Stack, F5 resume to normal completion | **PASS** |
| Pylance resolution | `import scrimp` no longer flagged after container `python.analysis.extraPaths` | **PASS** |
| Python interpreter GUI | `/opt/conda/envs/scrimp/bin/python` recognized after opting out of `Python Environments` manager | **PASS (user-confirmed)** |

### Required Dev Container settings

Use a **literal complete** image reference (the old
`${localEnv:SCRIMP_DOCKER_IMAGE:scrimp-mpi:v1.3.0-rc1}` expression was truncated
by VS Code and attempted to pull `scrimp-mpi:latest`). In
`project-template/.devcontainer/devcontainer.json`:

```json
{
  "image": "scrimp-mpi:v1.3.0-rc1",
  "customizations": {
    "vscode": {
      "settings": {
        "python.defaultInterpreterPath": "/opt/conda/envs/scrimp/bin/python",
        "python.useEnvironmentsExtension": false,
        "python.analysis.extraPaths": ["/workspace/scrimp"]
      }
    }
  }
}
```

This is a **partial excerpt**, not the full starter configuration.

The Python extension logged the correct active interpreter but
`ms-python.vscode-python-envs` reported `No Python environment is set for this
resource` and repeatedly offered to install an unrelated `uv` Python.
`python-envs.globalSearchPaths` did not solve this; editing JSON with a missing
comma also initially invalidated one local test. The developer
**confirmed the warnings ceased when the existing project actually contained**
`"python.useEnvironmentsExtension": false`. This setting affects the editor
environment manager **only**: do not install another Python, suppress missing
imports or alter the qualified numerical image. Pylance's independent false
`import scrimp` warning was fixed by `python.analysis.extraPaths`, with no
diagnostic suppression.

**Important:** an existing user project has its own
`.devcontainer/devcontainer.json` and is not automatically updated by a
`git pull` in the SCRIMP library repository. The starter template contains
the tested settings for newly created projects.

## Reproduction commands (inside VS Code Dev Container)

```bash
python -c "import sys, scrimp; print(sys.executable); print('SCRIMP OK')"
python simulations/first_simulation.py
mpiexec -n 2 python simulations/mpi_check.py
```

In VS Code: open `simulations/first_simulation.py`, set a breakpoint at
`system = S.DPHS("real")`, select `SCRIMP: first simulation`, press F5,
step with F10, inspect `system` and `q`, and press F5 to resume to completion.

## Acceptance boundaries: NOT release-qualified yet

- The tested image was built from `96961d64`. The subsequent source-level
  `DPHS._sync_lumped_trial_state` ComplexWarning correction **is not in that
  image**. Rebuild the candidate using the new branch HEAD and repeat relevant
  smoke checks without invalidating the native dependency cache.
- The PETSc unused `-info` option warning remains to be triaged.
- A **fresh-user / clean-room** onboarding run with a generated project,
  after fresh image pull, is recommended prior to public distribution.
- No qualification yet for Windows/WSL2, macOS or Apple Silicon image builds.
- The 2-rank starter only validates MPI startup; it is **not** MPI PDE physics
  validation, power-balance proof, scaling evidence or multiprocess debugging.
- Full 1/2/4/6-rank distributed PDE equivalence, discretized port-Hamiltonian
  power balances, Maxwell 3D and long-time memory/scaling are distinct roadmap
  acceptance gates.
- The image has **not** been pushed to GHCR, promoted as `latest`, or
  substituted for the developer's qualified `scrimp-mpi:local`.

## Next workstream

1. Rebuild `scrimp-mpi:v1.3.0-rc1` from the updated GitHub branch, measure
   warm-cache behavior and verify the `ComplexWarning` fix.
2. Run `scripts/smoke-v130-docker-candidate.sh` again and retain full logs
   and exact image revision.
3. Follow the remaining scientific validation matrix and performance plan in
   [V130_RELEASE_PLAN.md](V130_RELEASE_PLAN.md).
4. Only after the remaining scientific, platform and packaging gates: publish
   the tested versioned GHCR image and matching independent user-project ZIP.

**Conclusion:** the Ubuntu `SCRIMP + Docker + VS Code` developer workflow
is accepted; **the v1.3.0 release as a whole is not yet qualified**.
