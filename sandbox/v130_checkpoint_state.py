# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See LICENSE for license information.

"""#29: versioned, fail-closed accepted-state checkpoint *sandbox pilot*.

Pure NumPy/JSON I/O: no PETSc/MUMPS objects, pickle or factor data.
This is a COLD state restart, NOT an integrator-history checkpoint.
The pilot deliberately rejects a changed FE mesh, DOF ordering,
weak form, parameter, control or original physical time origin.

No general production DPHS checkpoint API is implied by this module.
"""

import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np

SCHEMA = "scrimp-v130-accepted-dae-state-checkpoint-v1"
RESTART_KIND = "cold_state_only_not_ts_history"
MAX_METADATA_BYTES = 4_000_000


def canonical(data):
    return json.dumps(
        data, sort_keys=True, ensure_ascii=True,
        separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def sha256(payload):
    return hashlib.sha256(payload).hexdigest()


def _digest_array(array):
    """Hash exact shape/type/ordering/bytes, no tolerance or sorting."""
    a = np.asarray(array)
    if a.dtype.kind not in "biufc":
        raise ValueError("Mesh/FEM fingerprint refuses nonnumeric GetFEM data")
    v = np.ascontiguousarray(a)
    return sha256(canonical({
        "shape": list(v.shape),
        "dtype": v.dtype.str,
    }) + v.tobytes(order="C"))


def model_provenance(system):
    """Fingerprint actual assembled GetFEM geometry and source weak forms.

    Includes per-convex geometric coordinates (connectivity/order),
    region identifiers, and ordered basic FEM DOF nodes. If a required
    GetFEM 5.5 API is absent, fail closed instead of hashing only counts.
    Actual mesh file persistence is a separate future product feature.
    """
    if system.domain is None:
        raise ValueError("Cannot checkpoint a system without a GetFEM domain")
    mesh_data = []
    for mesh in system.domain.get_mesh():
        pts_from_cv = mesh.pts_from_cvid()
        if not isinstance(pts_from_cv, (tuple, list)) or len(pts_from_cv) != 2:
            raise ValueError("GetFEM pts_from_cvid did not expose convex connectivity")
        mesh_data.append({
            "points": _digest_array(mesh.pts()),
            "cvid": _digest_array(mesh.cvid()),
            "convex_points": _digest_array(pts_from_cv[0]),
            "convex_offsets": _digest_array(pts_from_cv[1]),
            "dimension": int(mesh.dim()),
        })
    ports = []
    for name, port in system.ports.items():
        fem = port.get_fem()
        params = []
        for parameter in port.get_parameters():
            params.append({
                "name": parameter.get_name(),
                "kind": parameter.get_kind(),
                "expression": parameter.get_expression(),
            })
        fem_dof_connectivity = None
        fem_dof_offsets = None
        fem_description = None
        if fem is not None:
            # Coordinates alone are insufficient: FE orientation and
            # per-convex DOF numbering may differ with the same nodes.
            native_connectivity = fem.basic_dof_from_cvid()
            if (not isinstance(native_connectivity, (tuple,list))
                    or len(native_connectivity) != 2):
                raise ValueError(
                    "GetFEM did not expose ordered per-convex FE DOFs"
                )
            fem_dof_connectivity = _digest_array(native_connectivity[0])
            fem_dof_offsets = _digest_array(native_connectivity[1])
            fem_description = fem.char()
        ports.append({
            "name": str(name), "flow": port.get_flow(),
            "effort": port.get_effort(), "kind": port.get_kind(),
            "region": port.get_region(), "mesh_id": port.get_mesh_id(),
            "algebraic": port.get_algebraic(),
            "substituted": port.get_substituted(),
            "fem_nodes": (
                _digest_array(fem.basic_dof_nodes()) if fem is not None else None
            ),
            "fem_descriptor": fem_description,
            "fem_convex_dofs": fem_dof_connectivity,
            "fem_convex_dof_offsets": fem_dof_offsets,
            "fem_nbdof": int(fem.nbdof()) if fem is not None else None,
            "parameters": params,
        })
    bricks = []
    for name, brick in system.bricks.items():
        bricks.append({
            "name": str(name), "expression": brick.get_form(),
            "position": brick.get_position(),
            "dt": brick.get_dt(), "linear": brick.get_linear(),
            "explicit": brick.get_explicit(),
            "regions": list(brick.get_regions()),
            "mesh_id": brick.get_mesh_id(),
        })
    variables = []
    # Order is essential; an FE permutation must invalidate the state.
    for port in system.ports.values():
        for candidate in (port.get_flow(), port.get_effort()):
            if any(item["name"] == candidate for item in variables):
                continue
            try:
                start, width = system.gf_model.interval_of_variable(candidate)
            except Exception:
                continue  # substituted effort: not a model unknown
            variables.append({
                "name": candidate, "start": int(start), "width": int(width),
            })
    variables.sort(key=lambda item: (item["start"], item["name"]))
    n = int(system.gf_model.nbdof())
    if not variables or n <= 0:
        raise ValueError("No full ordered GetFEM state layout to checkpoint")
    domain = {
        "name": system.domain.get_name(),
        "dimensions": list(system.domain.get_dim()),
        "subdomains": system.domain.get_subdomains(),
        "boundaries": system.domain.get_boundaries(),
    }
    detail = {
        "ndof": n, "domain": domain, "meshes": mesh_data,
        "variables": variables, "ports": ports, "bricks": bricks,
    }
    return {
        "fingerprint_sha256": sha256(canonical(detail)),
        "detail": detail,
    }


def _validate_state(z):
    a = np.asarray(z)
    if a.ndim != 1 or a.size <= 0 or a.dtype.kind != "f":
        raise ValueError("Full accepted GetFEM state must be a real 1-D vector")
    if not np.all(np.isfinite(a)):
        raise ValueError("Checkpoint contains NaN or infinite FE coefficients")
    return np.asarray(a, dtype="<f8").copy()


def write_checkpoint(path, *, z, t, provenance, ts_type, accepted, extra=None):
    """Atomically save a previously accepted complete state on rank 0.

    Callers must certify that z is the actual accepted PETSc TS state;
    the function refuses unconverged or incomplete source runs.
    """
    if accepted is not True:
        raise ValueError("Refusing checkpoint of non-accepted PETSc TS state")
    arr = _validate_state(z)
    t = float(t)
    if not math.isfinite(t) or t < 0:
        raise ValueError("Invalid absolute physical checkpoint time")
    if ts_type != "cn":
        raise ValueError("Pilot supports only PETSc TSCN cold state checkpoints")
    if not isinstance(provenance, dict) or not isinstance(
        provenance.get("fingerprint_sha256"), str
    ):
        raise ValueError("Missing complete original GetFEM model fingerprint")
    if int(provenance.get("detail", {}).get("ndof", -1)) != arr.size:
        raise ValueError("Checkpoint full FE state and model DOF count mismatch")
    if sha256(canonical(provenance["detail"])) != provenance["fingerprint_sha256"]:
        raise ValueError("Corrupt original GetFEM model provenance")
    if extra is None:
        extra = {}
    if not isinstance(extra, dict):
        raise ValueError("Checkpoint extra metadata must be a JSON object")
    header = {
        "schema": SCHEMA, "kind": RESTART_KIND,
        "time_absolute": t, "scheme": "cn",
        "state_dtype": "<f8", "state_size": int(arr.size),
        "state_sha256": sha256(arr.tobytes(order="C")),
        "model": provenance, "extra": extra,
    }
    encoded = canonical(header)
    if len(encoded) > MAX_METADATA_BYTES:
        raise ValueError("Checkpoint metadata exceeds fixed schema size cap")
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, tempname = tempfile.mkstemp(
        prefix="." + destination.name + ".", suffix=".pending",
        dir=str(destination.parent),
    )
    try:
        with os.fdopen(fd, "wb") as handle:
            np.savez_compressed(
                handle, state=arr,
                manifest=np.frombuffer(encoded, dtype=np.uint8),
            )
            handle.flush()
            os.fsync(handle.fileno())
        # Atomic hard-link commit refuses overwriting previously accepted
        # scientific evidence. Source and destination share one filesystem.
        os.link(tempname, destination)
        os.unlink(tempname)
    finally:
        if os.path.exists(tempname):
            os.unlink(tempname)
    return header


def read_checkpoint(path, *, current_provenance, ts_type="cn"):
    """Validate full source and current model before returning any state."""
    with np.load(path, allow_pickle=False) as data:
        if set(data.files) != {"state", "manifest"}:
            raise ValueError("Incomplete or unknown checkpoint payload")
        raw = np.asarray(data["manifest"])
        if raw.ndim != 1 or raw.dtype != np.dtype("uint8"):
            raise ValueError("Invalid checkpoint metadata encoding")
        if raw.size > MAX_METADATA_BYTES:
            raise ValueError("Checkpoint manifest exceeds size cap")
        header = json.loads(raw.tobytes().decode("utf-8"))
        state = np.asarray(data["state"])
    if not isinstance(header, dict) or header.get("schema") != SCHEMA:
        raise ValueError("Unsupported checkpoint schema")
    if header.get("kind") != RESTART_KIND or header.get("scheme") != "cn":
        raise ValueError("Unsupported integrator history or checkpoint semantics")
    if ts_type != "cn":
        raise ValueError("Checkpoint cold restart requires original PETSc CN")
    if header.get("state_dtype") != "<f8" or state.dtype != np.dtype("<f8"):
        raise ValueError("Unexpected checkpoint scalar precision/type")
    arr = _validate_state(state)
    if header.get("state_size") != arr.size or (
        sha256(arr.tobytes(order="C")) != header.get("state_sha256")
    ):
        raise ValueError("Damaged checkpoint FE state hash or dimensions")
    t = header.get("time_absolute")
    if type(t) not in (float, int) or not math.isfinite(t) or t < 0:
        raise ValueError("Damaged physical checkpoint time")
    checkpoint_provenance = header.get("model", {})
    if not isinstance(checkpoint_provenance, dict):
        raise ValueError("Damaged model provenance structure")
    current = current_provenance
    if (not isinstance(current, dict)
            or current.get("fingerprint_sha256") !=
            checkpoint_provenance.get("fingerprint_sha256")
            or sha256(canonical(current.get("detail"))) !=
            current.get("fingerprint_sha256")
            or sha256(canonical(checkpoint_provenance.get("detail"))) !=
            checkpoint_provenance.get("fingerprint_sha256")):
        raise ValueError(
            "Checkpoint mesh, FE numbering, weak forms, physical parameters "
            "or boundary control differs from the reconstructed model"
        )
    if arr.size != int(current["detail"]["ndof"]):
        raise ValueError("Loaded FE vector size differs from current model")
    return float(t), arr, header


def collective_read(path, *, mpi, current_provenance, ts_type="cn"):
    """Rank-0 disk read and cross-rank collective fail-closed validation."""
    local_fp = current_provenance.get("fingerprint_sha256")
    fps = mpi.allgather(local_fp)
    if any(fp != fps[0] for fp in fps):
        raise ValueError("GetFEM model fingerprint disagrees across MPI ranks")
    payload = None
    if mpi.rank == 0:
        try:
            payload = (
                True, read_checkpoint(
                    path, current_provenance=current_provenance,
                    ts_type=ts_type,
                )
            )
        except Exception as err:
            payload = (False, f"{type(err).__name__}: {err}")
    ok, value = mpi.bcast(payload, root=0)
    if not ok:
        raise ValueError("Collective checkpoint import rejected: " + value)
    return value


def collective_write(path, *, mpi, z, t, provenance, ts_type="cn",
                     accepted, extra=None):
    """Persist only if *every* rank sees same accepted state and model."""
    # Ordinary per-rank validation failures must be collective:
    # one rank must not throw before its peers enter allgather.
    arr = None
    try:
        arr = _validate_state(z)
        evidence = (
            True, sha256(arr.tobytes(order="C")),
            provenance.get("fingerprint_sha256"), float(t), bool(accepted),
        )
    except Exception as error:
        evidence = (False, f"{type(error).__name__}: {error}")
    hashes = mpi.allgather(evidence)
    if any(item[0] is not True for item in hashes):
        raise ValueError("MPI checkpoint source rejected: "+repr(hashes))
    if any(item != hashes[0] for item in hashes):
        raise ValueError("MPI inconsistent accepted state/time/model for checkpoint")
    result = None
    if mpi.rank == 0:
        try:
            result = (
                True, write_checkpoint(
                    path, z=arr, t=t, provenance=provenance,
                    ts_type=ts_type, accepted=accepted, extra=extra,
                )
            )
        except Exception as err:
            result = (False, f"{type(err).__name__}: {err}")
    ok, value = mpi.bcast(result, root=0)
    if not ok:
        raise ValueError("Collective checkpoint export rejected: " + value)
    return value
