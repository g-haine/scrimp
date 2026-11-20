#!/usr/bin/env python3
# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Source-v5.5 GetFEM finite-element *constructor* audit.

Every probe runs in a separate Python process: GetFEM C++ assertions may
terminate the interpreter. Passing this audit means that a descriptor was
constructed with the expected number of local DoFs. It does *not* certify
Piola transforms, facet continuity, unisolvance, convergence or MPI behavior.

Examples:
    # Run INSIDE the scrimp-mpi:local GetFEM 5.5 container, NOT host conda.
    python sandbox/audit_getfem55_element_catalog.py --max-order 2
    python sandbox/audit_getfem55_element_catalog.py --max-order 5 --timeout 120
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

_PREFIX = "GETFEM_AUDIT_JSON:"


def expected_dofs(family: str, dim: int, order: int) -> int:
    if family == "RT":
        if dim == 2:
            return (order + 1) * (order + 3)
        return (order + 1) * (order + 2) * (order + 4) // 2
    if family == "RT0Q":
        return 2 * dim
    if family == "BDM":
        if dim == 2:
            return (order + 1) * (order + 2)
        return (order + 1) * (order + 2) * (order + 3) // 2
    if family == "N1":
        if dim == 2:
            return (order + 1) * (order + 3)
        return (order + 1) * (order + 3) * (order + 4) // 2
    raise ValueError(f"Unknown family: {family}")


def child_probe(name: str) -> int:
    """Only imported in the subprocess; avoid taking the parent down."""
    try:
        import getfem

        fem = getfem.Fem(name)
        result = {
            "descriptor": name,
            "nbdof": int(fem.nbdof()),
            "dim": int(fem.dim()),
            "target_dim": int(fem.target_dim()),
            "getfem_module": str(getfem.__file__),
        }
        print(_PREFIX + json.dumps(result, sort_keys=True), flush=True)
        return 0
    except BaseException as exc:
        print(
            _PREFIX
            + json.dumps(
                {"descriptor": name, "error": f"{type(exc).__name__}: {exc}"},
                sort_keys=True,
            ),
            flush=True,
        )
        return 2


def catalog(max_order: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for dim in (2, 3):
        for order in range(max_order + 1):
            rows.append(
                {
                    "family": "RT",
                    "dim": dim,
                    "order": order,
                    "descriptor": f"FEM_RTK({dim},{order})",
                    "expected_dofs": expected_dofs("RT", dim, order),
                }
            )
        for order in range(1, max_order + 1):
            rows.append(
                {
                    "family": "BDM",
                    "dim": dim,
                    "order": order,
                    "descriptor": f"FEM_BDMK({dim},{order})",
                    "expected_dofs": expected_dofs("BDM", dim, order),
                }
            )
        rows.append(
            {
                "family": "RT0Q",
                "dim": dim,
                "order": 0,
                "descriptor": f"FEM_RT0Q({dim})",
                "expected_dofs": expected_dofs("RT0Q", dim, 0),
            }
        )
        rows.append(
            {
                "family": "N1",
                "dim": dim,
                "order": 0,
                "descriptor": f"FEM_NEDELEC({dim})",
                "expected_dofs": expected_dofs("N1", dim, 0),
            }
        )
    return rows


def run_one(row: dict[str, Any], timeout: float) -> dict[str, Any]:
    cmd = [sys.executable, str(Path(__file__).resolve()), "--probe", row["descriptor"]]
    try:
        p = subprocess.run(cmd, text=True, capture_output=True, timeout=timeout, check=False)
    except subprocess.TimeoutExpired as exc:
        return {
            **row, "status": "timeout", "timeout_seconds": timeout,
            "stdout": (exc.stdout or b"").decode(errors="replace")[-1500:]
            if isinstance(exc.stdout, bytes) else (exc.stdout or "")[-1500:],
        }
    payload: dict[str, Any] = {}
    for line in p.stdout.splitlines():
        if line.startswith(_PREFIX):
            try:
                payload = json.loads(line[len(_PREFIX):])
            except json.JSONDecodeError:
                payload = {"invalid_json": line}
    observed = payload.get("nbdof")
    status = (
        "constructed_matching_dimension"
        if p.returncode == 0 and observed == row["expected_dofs"]
        and payload.get("dim") == row["dim"]
        and payload.get("target_dim") == row["dim"]
        else "dimension_mismatch" if p.returncode == 0
        else "constructor_error"
    )
    return {
        **row, **payload, "status": status, "returncode": p.returncode,
        "stdout_tail": p.stdout[-1200:],
        "stderr_tail": p.stderr[-2400:],
    }


def runtime_provenance(*, allow_non_docker: bool) -> dict[str, Any]:
    """Fail closed before probing a host-installed or serial GetFEM.

    This proves that the Python binding is loaded from the expected Docker
    environment and reports GetFEM MPI level; the image's source SHA and
    tarball checksum must still be recorded externally at build time.
    """
    try:
        import getfem
    except (ImportError, OSError) as exc:
        raise RuntimeError(f"Cannot load GetFEM in active Python: {exc}") from exc
    module = str(Path(getfem.__file__).resolve())
    mpi_level = int(getfem.util_mpi_parallelism_level())
    container_marker = Path("/.dockerenv").exists()
    expected_prefix = "/opt/conda/envs/scrimp/"
    provenance = {
        "python_executable": sys.executable,
        "python_prefix": sys.prefix,
        "getfem_module": module,
        "getfem_mpi_parallelism_level": mpi_level,
        "docker_marker_present": container_marker,
        "active_conda_prefix": os.environ.get("CONDA_PREFIX"),
        "source_commit_verified": False,
        "note": "Image source identity requires external build provenance.",
    }
    if not allow_non_docker:
        issues = []
        if not container_marker:
            issues.append("not running inside a Docker container")
        if not str(Path(sys.executable).resolve()).startswith(expected_prefix):
            issues.append(f"Python executable is not in {expected_prefix}")
        if not module.startswith(expected_prefix):
            issues.append(f"GetFEM module is not in {expected_prefix}")
        if mpi_level != 2:
            issues.append(f"GetFEM MPI level is {mpi_level}, expected 2")
        if issues:
            raise RuntimeError(
                "REFUSING non-qualified GetFEM constructor audit: "
                + "; ".join(issues)
                + ". Use the SCRIMP GetFEM-5.5 Docker image. "
                "For deliberate host diagnostics only, pass --allow-non-docker "
                "(results will be flagged as NON-QUALIFYING)."
            )
    provenance["environment_qualified_for_docker_probe"] = not allow_non_docker
    return provenance


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-order", type=int, default=2,
                        help="Highest RT/BDM polynomial order to probe (default 2)")
    parser.add_argument("--timeout", type=float, default=90.0,
                        help="Per-element subprocess timeout, in seconds")
    parser.add_argument("--output", type=Path,
                        default=Path("sandbox/outputs/getfem55_constructor_audit.json"),
                        help="JSON report destination; use '-' for JSON on stdout "
                             "(status/progress goes to stderr)")
    parser.add_argument("--only", action="append", metavar="DESCRIPTOR",
                        help="Probe one exact descriptor; repeat for several "
                             "(for example --only 'FEM_BDMK(3,1)')")
    parser.add_argument("--strict", action="store_true",
                        help="Exit 1 if any constructor did not match its expected dimensions")
    parser.add_argument("--allow-non-docker", action="store_true",
                        help="Diagnostic only: override the mandatory Docker/GetFEM MPI-2 guard")
    parser.add_argument("--probe", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.probe is not None:
        return child_probe(args.probe)
    json_to_stdout = args.output == Path("-")
    progress_stream = sys.stderr if json_to_stdout else sys.stdout

    def progress(*items: object, **kwargs: Any) -> None:
        print(*items, file=progress_stream, **kwargs)

    if not (0 <= args.max_order <= 5):
        parser.error("--max-order must be between 0 and 5")
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    try:
        provenance = runtime_provenance(allow_non_docker=args.allow_non_docker)
    except RuntimeError as exc:
        parser.error(str(exc))
    progress("Runtime provenance:", json.dumps(provenance, sort_keys=True), flush=True)
    rows = catalog(args.max_order)
    if args.only:
        available = {row["descriptor"] for row in rows}
        missing = set(args.only) - available
        if missing:
            parser.error(f"Requested descriptor(s) outside --max-order catalog: {sorted(missing)}")
        chosen = set(args.only)
        rows = [row for row in rows if row["descriptor"] in chosen]
    results = []
    for row in rows:
        progress(f"{row['descriptor']:<22} expected={row['expected_dofs']:<4}", end=" ", flush=True)
        result = run_one(row, args.timeout)
        results.append(result)
        observed = result.get("nbdof", "-")
        progress(f"observed={str(observed):<4} {result['status']}", flush=True)
        if result["status"] != "constructed_matching_dimension":
            why = result.get("error") or result.get("stderr_tail") or "no diagnostic"
            progress("  detail:", str(why).strip()[:1200], flush=True)
    out = {
        "audit_type": "constructor_and_local_dimension_only",
        "reference_source_expected": "GetFEM v5.5 tag c242f7fe0989e9a90846c63eaccc6b495d774724",
        "runtime_provenance": provenance,
        "qualifying_environment": not args.allow_non_docker,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python_executable": sys.executable,
        "max_order": args.max_order,
        "only_descriptors": args.only,
        "results": results,
    }
    json_report = json.dumps(out, indent=2, ensure_ascii=False) + "\n"
    if json_to_stdout:
        sys.stdout.write(json_report)
        sys.stdout.flush()
    else:
        try:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json_report, encoding="utf-8")
        except OSError as exc:
            print(
                f"WARNING: unable to write JSON to {args.output}: {exc}. "
                "Emitting full report on stdout between markers instead; "
                "rerun with --output - and redirect stdout on the host.",
                file=sys.stderr,
                flush=True,
            )
            print("GETFEM_AUDIT_FALLBACK_JSON_BEGIN", flush=True)
            sys.stdout.write(json_report)
            print("GETFEM_AUDIT_FALLBACK_JSON_END", flush=True)
    bad = sum(row["status"] != "constructed_matching_dimension" for row in results)
    progress(f"\nConstructor audit: {len(results) - bad}/{len(results)} constructed "
             f"with expected dimensions; {bad} errors/mismatches/timeouts.")
    progress(f"JSON report: {'stdout' if json_to_stdout else args.output}")
    progress("WARNING: MPI level 2 is checked, but this is NOT an MPI-parity, "
             "Piola or H(div)/H(curl) validation.")
    return int(args.strict and bad > 0)


if __name__ == "__main__":
    sys.exit(main())
