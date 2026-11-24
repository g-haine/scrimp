# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Validate rank-0-only MPI ParaView field data, geometry and PVD series.

This is a **host-side audit** of real native GetFEM VTK output. The
only verified scalar formats are the original p and T fields. The
GetFEM XML writer can emit either a standard 4-byte binary length
prefix or a float32 point count for the scalar PointData array.

Do not confuse validity of one global rank-0 piece with a generic
domain-decomposed .pvtu exporter. This accepts the Heat-Wave case only
after full MPI2 FEM state and geometric coverage were qualified.
"""

import base64
import hashlib
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np


def _bytes_array(data_array):
    if data_array is None or data_array.get("format") != "binary":
        raise RuntimeError("Expected actual GetFEM binary VTK DataArray")
    return base64.b64decode((data_array.text or "").strip(), validate=True)


def _scalar_field(blob, field):
    piece = ET.fromstring(blob).find(".//Piece")
    if piece is None:
        raise RuntimeError("Missing VTK UnstructuredGrid Piece")
    points = int(piece.attrib["NumberOfPoints"])
    cells = int(piece.attrib["NumberOfCells"])
    if points < 1 or cells < 1:
        raise RuntimeError("Missing complete mesh points or cells")

    arrays = {}
    for group, names in (
        ("Points", ("Points",)),
        ("Cells", ("connectivity", "offsets", "types")),
    ):
        element = piece.find(group)
        if element is None:
            raise RuntimeError(f"VTK mesh lacks {group}")
        for name in names:
            data = next(
                (a for a in element.findall("DataArray")
                 if a.attrib.get("Name") == name), None
            )
            arrays[group + "/" + name] = _bytes_array(data)

    values_node = piece.find("PointData")
    if values_node is None:
        raise RuntimeError("VTK missing PointData")
    scalar = next(
        (a for a in values_node.findall("DataArray")
         if a.attrib.get("Name") == field), None
    )
    if scalar is None or scalar.attrib.get("type") != "Float32":
        raise RuntimeError(f"Missing scalar Float32 field {field}")
    raw = _bytes_array(scalar)
    if len(raw) != 4*(points+1):
        raise RuntimeError("Unexpected scalar VTK payload size")
    payload = np.frombuffer(raw, dtype="<f4")
    standard_length = int(np.frombuffer(raw[:4], dtype="<u4")[0])
    getfem_point_count = float(payload[0])
    if standard_length != points*4 and getfem_point_count != float(points):
        raise RuntimeError("Invalid GetFEM scalar array binary prefix")
    values = payload[1:].copy()
    if not np.all(np.isfinite(values)):
        raise RuntimeError("Nonfinite VTK field values")
    metadata = {
        "points": points,
        "cells": cells,
        "field_min": float(np.min(values)),
        "field_max": float(np.max(values)),
        "field_l2": float(np.linalg.norm(values.astype(np.float64))),
        "mesh_arrays_sha256": hashlib.sha256(
            b"".join(arrays[key] for key in sorted(arrays))
        ).hexdigest(),
        "field_sha256": hashlib.sha256(values.tobytes()).hexdigest(),
    }
    return metadata, values


def audit_root_only_pv(case_dir):
    """Validate serial and two-rank shared-PV output under one heat case.

    case_dir: directory containing heatX_mpi1/ and heatX_mpi2/.
    All expected output paths are user-supplied artifact data.
    """
    case_dir = Path(case_dir)
    results = {}
    for heat in (1, 2):
        serial_dir = case_dir / f"heat{heat}_mpi1" / "pv"
        parallel_dir = case_dir / f"heat{heat}_mpi2" / "pv"
        for field in ("p", "T"):
            series = {}
            for label, folder in (
                ("MPI1", serial_dir / field),
                ("MPI2", parallel_dir / field),
            ):
                collection = ET.parse(folder / f"{field}.pvd")
                datasets = collection.findall(".//DataSet")
                if len(datasets) != 5:
                    raise RuntimeError(f"Expected five {label} PVD steps for {field}")
                snapshots = []
                for step, item in enumerate(datasets):
                    t = float(item.attrib["timestep"])
                    if not np.isclose(t, step*0.001, rtol=0, atol=1e-12):
                        raise RuntimeError("Incorrect PVD physical timestep")
                    name = item.attrib["file"]
                    if Path(name).name != name:
                        raise RuntimeError("Unsafe PVD file path")
                    content = (folder / name).read_bytes()
                    metadata, values = _scalar_field(content, field)
                    metadata["time"] = t
                    metadata["bytes"] = len(content)
                    metadata["sha256"] = hashlib.sha256(content).hexdigest()
                    snapshots.append((metadata, values))
                series[label] = snapshots
            comparisons = []
            for (ref_meta, ref), (cur_meta, cur) in zip(
                series["MPI1"], series["MPI2"]
            ):
                same_geometry = (
                    ref_meta["mesh_arrays_sha256"] ==
                    cur_meta["mesh_arrays_sha256"]
                    and ref_meta["points"] == cur_meta["points"]
                    and ref_meta["cells"] == cur_meta["cells"]
                )
                max_difference = float(np.max(np.abs(ref-cur)))
                field_equal = bool(
                    np.allclose(ref, cur, rtol=2e-6, atol=2e-7)
                )
                if not same_geometry or not field_equal:
                    raise RuntimeError(
                        f"Root-only MPI2 VTK mismatch: heat={heat}, "
                        f"field={field}, t={ref_meta['time']}: "
                        f"geometry={same_geometry}, data={max_difference}"
                    )
                comparisons.append({
                    "t": ref_meta["time"],
                    "points": cur_meta["points"],
                    "cells": cur_meta["cells"],
                    "max_abs_field_error": max_difference,
                    "same_mesh_binary_data": same_geometry,
                    "field_within_float32_tolerance": field_equal,
                    "MPI1_field_SHA256": ref_meta["field_sha256"],
                    "MPI2_field_SHA256": cur_meta["field_sha256"],
                })
            results[f"heat{heat}_{field}"] = comparisons
        # Rank-1 must neither create an accidental shadow PVD/VTU
        # nor write into a supposedly exclusive root-only directory.
        rank_one = case_dir / f"heat{heat}_mpi2" / "rank1"
        if rank_one.exists():
            raise RuntimeError(f"Unexpected rank1 ParaView output: {rank_one}")
    return results
