# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Host-side regression for the real MPI CN comparison's inline Python.

Exercises the *actual* inline comparator extracted from its shell
driver, rather than a rewritten imitation, with fake completed
per-rank reports. Does not claim to run or replace PETSc/GetFEM.
This test would have caught the historical uninitialized host_error
in the SUCCESS path, and checks archiving even when a comparison
fails or a native MPI rank exited nonzero.
"""

import base64
import contextlib
import io
import json
import os
import struct
from pathlib import Path
from types import SimpleNamespace
from xml.etree import ElementTree as ET
import sys
import tempfile
import traceback
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "run-v130-auto-interface-mpi-cn.sh"
START = 'python3 - "$out" <<\'PY\'\n'
STOP = '\nPY\n'


def inline_comparator():
    source = DRIVER.read_text(encoding="utf-8")
    if source.count(START) != 1:
        raise AssertionError("Expected one native-driver inline Python comparator")
    tail = source.split(START, 1)[1]
    if STOP not in tail:
        raise AssertionError("Comparator heredoc terminator is missing")
    return tail.split(STOP, 1)[0] + "\n"


def prepare_fake_evidence(root):
    """Generate minimal reports and deterministic *simulated* trace files."""
    for heat in (1, 2):
        for size in (1, 2):
            label = f"heat{heat}_mpi{size}"
            folder = root / label
            folder.mkdir()
            (root / f"{label}.exitcode").write_text("0\n")
            for rank in range(size):
                record = {
                    "MPI_ranks": size,
                    "heat_region": heat,
                    "rank": rank,
                    "status": "PASS_RANK_MPI_CN_AND_PHYSICAL_GATES",
                    "rank_consensus": True,
                    "real_PETSc_CN_steps": 4,
                    "physical_gates": {"fake_physics": True},
                    "plot_and_ParaView_export_skipped_on_all_ranks": True,
                    "MPI_rank_isolated_real_ParaView_export_diagnostic": False,
                    "MPI_root_only_shared_ParaView_export_pilot": False,
                    "physical_GetFEM_assembly_executed_on_all_ranks": True,
                    "rank_fingerprints": [{"sha": "same"}] * size,
                    "automatic_interface": {"logical_region": 10},
                }
                (folder / f"rank{rank}.json").write_text(
                    json.dumps(record, indent=2) + "\n"
                )
                np.savez_compressed(
                    folder / f"rank{rank}.npz",
                    times=np.arange(5, dtype=float) / 1000,
                    states=np.arange(20, dtype=float).reshape(5, 4),
                    H=np.array([5., 4.9, 4.8, 4.7, 4.6]),
                    phys_CN_P_interface_heat=np.ones(5),
                    phys_CN_P_interface_wave=np.ones(5),
                )


def prepare_fake_pv_evidence(root):
    """Construct *synthetic* parsable PVD/VTU files for host logic tests.

    These files do not represent GetFEM meshes; only the native
    diagnostic can establish actual parallel VTK export correctness.
    """
    for heat in (1, 2):
        for size in (1, 2):
            label = f"heat{heat}_mpi{size}"
            for rank in range(size):
                base = root / label / f"rank{rank}"
                report_path = root / label / f"rank{rank}.json"
                report = json.loads(report_path.read_text())
                report["plot_and_ParaView_export_skipped_on_all_ranks"] = False
                report["MPI_rank_isolated_real_ParaView_export_diagnostic"] = True
                report_path.write_text(json.dumps(report) + "\n")
                for field in ("p", "T"):
                    folder = base / "pv" / field
                    folder.mkdir(parents=True)
                    pvd_root = ET.Element("VTKFile", type="Collection")
                    collection = ET.SubElement(pvd_root, "Collection")
                    for step in range(5):
                        name = f"{field}_{step:05d}.vtu"
                        ET.SubElement(
                            collection, "DataSet",
                            timestep=str(step*0.001), file=name,
                        )
                        vtk_root = ET.Element("VTKFile", type="UnstructuredGrid")
                        grid = ET.SubElement(vtk_root, "UnstructuredGrid")
                        ET.SubElement(
                            grid, "Piece",
                            NumberOfPoints="5", NumberOfCells="3",
                        )
                        (folder / name).write_bytes(
                            ET.tostring(vtk_root, encoding="utf-8")
                        )
                    (folder / f"{field}.pvd").write_bytes(
                        ET.tostring(pvd_root, encoding="utf-8")
                    )


def prepare_fake_root_only_pv_evidence(root):
    """Create tiny synthetic GetFEM-like scalar/mesh binary VTK arrays.

    This qualifies only HOST validation logic; the real native pilot
    must demonstrate full GetFEM mesh and accurate numerical values.
    """
    def array(parent, name, dtype, values, *, components=None, field=False):
        a = ET.SubElement(parent, "DataArray", Name=name, type=dtype,
                          format="binary")
        if components is not None:
            a.set("NumberOfComponents", str(components))
        values = np.asarray(values, dtype="<f4" if dtype == "Float32" else "<i4")
        prefix = (np.float32(len(values)).tobytes()
                  if field else struct.pack("<I", values.nbytes))
        a.text = base64.b64encode(prefix + values.tobytes()).decode("ascii")

    for heat in (1, 2):
        for size in (1, 2):
            label = f"heat{heat}_mpi{size}"
            # Existing report records were created by prepare_fake_evidence.
            for rank in range(size):
                report_path = root / label / f"rank{rank}.json"
                report = json.loads(report_path.read_text())
                report["plot_and_ParaView_export_skipped_on_all_ranks"] = False
                report["MPI_root_only_shared_ParaView_export_pilot"] = True
                report_path.write_text(json.dumps(report) + "\n")
            for field in ("p", "T"):
                folder = root / label / "pv" / field
                folder.mkdir(parents=True)
                collection_root = ET.Element("VTKFile", type="Collection")
                collection = ET.SubElement(collection_root, "Collection")
                for step in range(5):
                    filename = f"{field}_{step:05d}.vtu"
                    ET.SubElement(collection, "DataSet",
                                  timestep=str(step*0.001), file=filename)
                    vtk_root = ET.Element("VTKFile", type="UnstructuredGrid")
                    grid = ET.SubElement(vtk_root, "UnstructuredGrid")
                    piece = ET.SubElement(grid, "Piece",
                                          NumberOfPoints="5",
                                          NumberOfCells="3")
                    pts = ET.SubElement(piece, "Points")
                    array(pts, "Points", "Float32",
                          np.arange(15, dtype=np.float32), components=3)
                    cells = ET.SubElement(piece, "Cells")
                    array(cells, "connectivity", "Int32", np.arange(9))
                    array(cells, "offsets", "Int32", [3, 6, 9])
                    array(cells, "types", "Int32", [5, 5, 5])
                    pointdata = ET.SubElement(piece, "PointData")
                    data = np.full(5, step*.001+heat,
                                   dtype=np.float32)
                    array(pointdata, field, "Float32", data, field=True)
                    (folder / filename).write_bytes(ET.tostring(vtk_root))
                (folder / f"{field}.pvd").write_bytes(
                    ET.tostring(collection_root)
                )


def compare(root, *, pv_diagnostic=False, pv_root_only=False):
    """Execute the EXACT production inline comparator in this process.

    Starting another Python interpreter is unnecessary to test an inline
    script and has been observed to terminate with SIGKILL (-9) in the
    native conda environment. Running in-process avoids a second Python
    process (and an independent NumPy/BLAS import) while still exercising
    the actual comparator source, with its real sys.argv and exit code.
    All global state and captured streams are restored on exit.
    """
    stdout, stderr = io.StringIO(), io.StringIO()
    exit_code = 0
    namespace = {"__name__": "__main__", "__file__": str(DRIVER)}
    with patch.dict(os.environ, {
        "SCRIMP_V130_MPI_PV_DIAGNOSTIC": "1" if pv_diagnostic else "0",
        "SCRIMP_V130_MPI_PV_ROOT_ONLY": "1" if pv_root_only else "0",
    }):
        with patch.object(sys, "argv", [str(DRIVER), str(root)]):
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                try:
                    exec(compile(inline_comparator(), str(DRIVER), "exec"), namespace)
                except SystemExit as exc:
                    exit_code = int(exc.code) if isinstance(exc.code, int) else 1
                except Exception:
                    traceback.print_exc()
                    exit_code = 1
    return SimpleNamespace(
        returncode=exit_code, stdout=stdout.getvalue(),
        stderr=stderr.getvalue(),
    )


class HostComparatorRegression(unittest.TestCase):
    def test_actual_embedded_comparator_is_valid_python(self):
        compile(inline_comparator(), str(DRIVER), "exec")

    def test_success_writes_summary_and_archive_without_nameerror(self):
        with tempfile.TemporaryDirectory(prefix="v130-mpi-cn-success-") as d:
            root = Path(d) / "run"
            root.mkdir()
            prepare_fake_evidence(root)
            outcome = compare(root)
            self.assertEqual(outcome.returncode, 0, outcome.stderr + outcome.stdout)
            summary = json.loads((root / "comparison.json").read_text())
            self.assertEqual(
                summary["status"],
                "PASS_BOTH_AUTO_INTERFACE_ORIENTATIONS_MPI2_DAE_CN_ENERGY_PARITY",
            )
            self.assertIsNone(summary["host_comparison_failure"])
            self.assertEqual(set(summary["cases"]), {
                "heat1_wave2", "heat2_wave1"
            })
            self.assertTrue(Path(str(root) + ".zip").is_file())

    def test_real_trace_mismatch_fails_but_preserves_archive(self):
        with tempfile.TemporaryDirectory(prefix="v130-mpi-cn-diff-") as d:
            root = Path(d) / "run"
            root.mkdir()
            prepare_fake_evidence(root)
            path = root / "heat1_mpi2" / "rank0.npz"
            with np.load(path) as old:
                values = {key: np.asarray(old[key]).copy() for key in old.files}
            values["states"][1, 2] += .5
            np.savez_compressed(path, **values)
            outcome = compare(root)
            self.assertNotEqual(outcome.returncode, 0)
            summary = json.loads((root / "comparison.json").read_text())
            self.assertEqual(
                summary["status"],
                "FAIL_OR_UNQUALIFIED_AUTO_INTERFACE_MPI2_DAE_CN_PHYSICS",
            )
            self.assertEqual(
                summary["cases"]["heat1_wave2"]["status"],
                "FAIL_OR_UNQUALIFIED_MPI2_DAE_CN_PARITY",
            )
            self.assertIsNone(summary["host_comparison_failure"])
            self.assertTrue(Path(str(root) + ".zip").exists())

    def test_opt_in_pv_diagnostic_parses_all_synthetic_rank_files(self):
        with tempfile.TemporaryDirectory(prefix="v130-mpi-pv-success-") as d:
            root = Path(d) / "run"
            root.mkdir()
            prepare_fake_evidence(root)
            prepare_fake_pv_evidence(root)
            outcome = compare(root, pv_diagnostic=True)
            self.assertEqual(
                outcome.returncode, 0, outcome.stderr + outcome.stdout
            )
            report = json.loads((root / "comparison.json").read_text())
            self.assertEqual(
                report["status"],
                "DIAGNOSTIC_MPI_RANK_ISOLATED_PV_FILES_COMPLETE_CN_ENERGY_PASSED",
            )
            self.assertTrue(report["isolated_MPI_ParaView_diagnostic_executed"])
            self.assertFalse(report["shared_production_MPI_ParaView_qualified"])
            self.assertEqual(len(report["isolated_MPI_ParaView_evidence"]), 6)
            comparisons = report["rank_isolated_VTU_PVD_mesh_parity_by_field"]
            self.assertEqual(len(comparisons), 8)
            self.assertTrue(all(x["same_as_serial"] for x in comparisons.values()))
            self.assertTrue(Path(str(root) + ".zip").is_file())

    def test_pv_diagnostic_native_failure_keeps_specific_failure_status(self):
        with tempfile.TemporaryDirectory(prefix="v130-mpi-pv-fail-") as d:
            root = Path(d) / "run"
            root.mkdir()
            prepare_fake_evidence(root)
            prepare_fake_pv_evidence(root)
            (root / "heat2_mpi2.exitcode").write_text("1\n")
            outcome = compare(root, pv_diagnostic=True)
            self.assertNotEqual(outcome.returncode, 0)
            report = json.loads((root / "comparison.json").read_text())
            self.assertEqual(
                report["status"],
                "FAIL_OR_UNQUALIFIED_AUTO_INTERFACE_MPI2_DAE_CN_OR_PV_DIAGNOSTIC",
            )
            self.assertTrue(Path(str(root) + ".zip").is_file())

    def test_root_only_shared_vtu_host_mode_parses_real_scalar_layout(self):
        with tempfile.TemporaryDirectory(prefix="v130-mpi-pv-root-") as d:
            root = Path(d) / "run"
            root.mkdir()
            prepare_fake_evidence(root)
            prepare_fake_root_only_pv_evidence(root)
            outcome = compare(root, pv_root_only=True)
            self.assertEqual(outcome.returncode, 0, outcome.stderr+outcome.stdout)
            report = json.loads((root / "comparison.json").read_text())
            self.assertEqual(
                report["status"],
                "PASS_MPI2_ROOT_ONLY_SHARED_VTU_PVD_COMPLETE_CN_ENERGY",
            )
            self.assertTrue(report["shared_production_MPI_ParaView_qualified"])
            audit = report["shared_MPI_root_only_field_and_mesh_data_audit"]
            self.assertEqual(set(audit), {"heat1_p", "heat1_T", "heat2_p", "heat2_T"})
            self.assertTrue(all(
                all(row["field_within_float32_tolerance"]
                    and row["same_mesh_binary_data"] for row in rows)
                for rows in audit.values()
            ))
            self.assertTrue(Path(str(root)+".zip").is_file())

    def test_root_only_changed_field_fails_and_archives_evidence(self):
        with tempfile.TemporaryDirectory(prefix="v130-mpi-pv-root-bad-") as d:
            root = Path(d) / "run"
            root.mkdir()
            prepare_fake_evidence(root)
            prepare_fake_root_only_pv_evidence(root)
            path = root / "heat1_mpi2" / "pv" / "p" / "p_00004.vtu"
            xml = ET.fromstring(path.read_bytes())
            data = xml.find(".//PointData/DataArray")
            raw = base64.b64decode(data.text)
            values = np.frombuffer(raw, dtype="<f4").copy()
            values[-1] += 1.
            data.text = base64.b64encode(values.tobytes()).decode("ascii")
            path.write_bytes(ET.tostring(xml))
            outcome = compare(root, pv_root_only=True)
            self.assertNotEqual(outcome.returncode, 0)
            report = json.loads((root / "comparison.json").read_text())
            self.assertIn("Root-only MPI2 VTK mismatch",
                          report["host_comparison_failure"])
            self.assertTrue(Path(str(root)+".zip").is_file())

    def test_nonzero_native_rank_exit_fails_with_archive(self):
        with tempfile.TemporaryDirectory(prefix="v130-mpi-cn-native-fail-") as d:
            root = Path(d) / "run"
            root.mkdir()
            prepare_fake_evidence(root)
            (root / "heat2_mpi2.exitcode").write_text("1\n")
            outcome = compare(root)
            self.assertNotEqual(outcome.returncode, 0)
            summary = json.loads((root / "comparison.json").read_text())
            self.assertEqual(
                summary["cases"]["heat2_wave1"]["MPI2"]["exit_code"], 1,
            )
            self.assertIsNone(summary["host_comparison_failure"])
            self.assertTrue(Path(str(root) + ".zip").is_file())

    def test_malformed_rank_report_archived_as_host_error(self):
        with tempfile.TemporaryDirectory(prefix="v130-mpi-cn-bad-meta-") as d:
            root = Path(d) / "run"
            root.mkdir()
            prepare_fake_evidence(root)
            (root / "heat1_mpi2" / "rank1.json").write_text("not JSON\n")
            outcome = compare(root)
            self.assertNotEqual(outcome.returncode, 0)
            summary = json.loads((root / "comparison.json").read_text())
            self.assertIn("JSONDecodeError", summary["host_comparison_failure"])
            self.assertTrue(Path(str(root) + ".zip").is_file())


if __name__ == "__main__":
    unittest.main()
