# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.

"""#29: verify ONE native PRODUCTION DPHS checkpoint API phase.

This is not the final full 16-case comparator. It rejects a partial,
truncated or physically mismatched seed/continuous/restart archive,
including a checkpoint that differs from the *accepted* original TS
endpoint. Called both after Docker and before skipping an old case.
"""

import argparse
import json
from pathlib import Path
import re

import numpy as np

from compare_v130_checkpoint_public_api import segment
# Host-only audit must NOT import scrimp.__init__:
# scrimp initializes PETSc with sys.argv and treats --evidence/--mesh
# as unused PETSc options. Load this pure NumPy/JSON helper by exact file.
import importlib.util
ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "_scrimp_native_checkpoint_host_audit",
    ROOT / "scrimp" / "utils" / "dae_checkpoint.py",
)
_checkpoint_runtime = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_checkpoint_runtime)
read_checkpoint = _checkpoint_runtime.read_checkpoint
sha256 = _checkpoint_runtime.sha256

MESHES = (0.2, 0.1)
RANKS = (1, 2)
PHASES = ("seed", "continuous", "resume")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_accepted_checkpoint(root, mesh, ranks):
    label = f"h{str(mesh).replace('.', 'p')}"
    path = root / f"{label}_source_mpi{ranks}" / "accepted_state.npz"
    require(path.is_file(), f"Missing actual source accepted-state file {path}")
    # Explicit file ownership closes a corrupted NPZ even when
    # np.load raises before constructing an NpzFile.
    with path.open("rb") as handle:
        with np.load(handle, allow_pickle=False) as data:
            require(set(data.files) == {"state", "manifest"},
                    "Incomplete accepted-state checkpoint archive")
            raw = np.asarray(data["manifest"])
            require(raw.ndim == 1 and raw.dtype == np.dtype("uint8"),
                    "Invalid accepted-state manifest encoding")
            header = json.loads(raw.tobytes().decode("utf-8"))
    require(isinstance(header, dict) and isinstance(header.get("model"), dict),
            "Invalid checkpoint provenance manifest")
    time, values, verified = read_checkpoint(
        path, current_provenance=header["model"], ts_type="cn"
    )
    require(abs(time-.02) <= 1e-12,
            "Incorrect checkpoint absolute physical time")
    require(verified.get("extra", {}).get("source_mpi_ranks") == ranks
            and verified["extra"].get("nominal_mesh_size") == mesh
            and verified["extra"].get("physical_control_clock") == "absolute"
            and verified["extra"].get("physical_rotation_target") == .2,
            "Incorrect checkpoint original physical source/provenance")
    return values, verified


def verify_one(root, mesh, phase, ranks, source_ranks=None):
    root = Path(root)
    require(mesh in MESHES and phase in PHASES and ranks in RANKS,
            "Unsupported checkpoint native phase or mesh")
    if phase == "resume":
        require(source_ranks in RANKS,
                "Restart phase must specify checkpoint source MPI count")
    else:
        require(source_ranks is None,
                "Only restarted phases may identify source MPI")
    result = segment(root, mesh, phase, ranks, source_rank=source_ranks)
    r = result["report"]
    fingerprint = r.get("model_fingerprint_sha256", "")
    require(isinstance(fingerprint, str) and
            re.fullmatch(r"[a-f0-9]{64}", fingerprint) is not None,
            "Native original GetFEM model fingerprint missing")

    label = f"h{str(mesh).replace('.', 'p')}"
    case_name = (f"{label}_restart_from{source_ranks}_to{ranks}"
                 if phase == "resume" else f"{label}_{phase}_mpi{ranks}")
    log = root / f"{case_name}.log"
    require(log.is_file(), f"Native execution log missing for {case_name}")
    text = log.read_text(encoding="utf-8", errors="replace")
    require(
        f"CHECKPOINT PUBLIC API PHASE={phase} mesh={mesh:g} MPI={ranks} " in text,
        f"Native physical phase completed but expected success marker absent: {case_name}"
    )
    if phase == "seed":
        state, header = load_accepted_checkpoint(root, mesh, ranks)
        require(np.array_equal(state, result["states"][-1]),
                "Accepted FE checkpoint differs from original PETSc final state")
        require(header["model"]["fingerprint_sha256"] == fingerprint,
                "Checkpoint model fingerprint differs from seed original GetFEM")
        require(r.get("complete_final_state_sha256") ==
                sha256(result["states"][-1].tobytes()),
                "Original seed final FE vector hash differs from report")
    elif phase == "resume":
        state, header = load_accepted_checkpoint(root, mesh, source_ranks)
        require(np.array_equal(result["states"][0], state),
                "Native cold restart changed original accepted FE coefficients")
        require(header["model"]["fingerprint_sha256"] == fingerprint,
                "Restart reconstructed a different GetFEM mesh or physical model")
        restored = r["initialization_report"]
        require(restored.get("status")==
                "COMPLETE_DAE_STATE_COLD_RESTORED_IN_DPHS"
                and restored.get(
                    "all_GetFEM_and_GlobalVec_coefficients_restored_bitwise"
                ) is True
                and restored.get("model_sha256") ==
                    header["model"]["fingerprint_sha256"],
                "Native production DPHS publication differs from accepted state")
    else:
        require(r.get("complete_final_state_sha256") ==
                sha256(result["states"][-1].tobytes()),
                "Original continuous final FE state hash differs from report")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--evidence", type=Path, required=True)
    p.add_argument("--mesh", type=float, choices=MESHES, required=True)
    p.add_argument("--phase", choices=PHASES, required=True)
    p.add_argument("--ranks", type=int, choices=RANKS, required=True)
    p.add_argument("--source-ranks", type=int, choices=RANKS)
    args = p.parse_args()
    verify_one(
        args.evidence, args.mesh, args.phase, args.ranks,
        source_ranks=args.source_ranks,
    )
    print(f"PASS: verified native checkpoint phase {args.phase} "
          f"mesh={args.mesh:g} MPI={args.ranks}", flush=True)


if __name__ == "__main__":
    main()
