#!/usr/bin/env python3
"""Source-v5.5 GetFEM finite-element *constructor* audit.

Every probe runs in a separate Python process: GetFEM C++ assertions may
terminate the interpreter. Passing this audit means that a descriptor was
constructed with the expected number of local DoFs. It does *not* certify
Piola transforms, facet continuity, unisolvance, convergence or MPI behavior.

Examples:
    python sandbox/audit_getfem55_element_catalog.py --max-order 2
    python sandbox/audit_getfem55_element_catalog.py --max-order 5 --timeout 120
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-order", type=int, default=2,
                        help="Highest RT/BDM polynomial order to probe (default 2)")
    parser.add_argument("--timeout", type=float, default=90.0,
                        help="Per-element subprocess timeout, in seconds")
    parser.add_argument("--output", type=Path,
                        default=Path("sandbox/outputs/getfem55_constructor_audit.json"))
    parser.add_argument("--strict", action="store_true",
                        help="Exit 1 if any constructor did not match its expected dimensions")
    parser.add_argument("--probe", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.probe is not None:
        return child_probe(args.probe)
    if not (0 <= args.max_order <= 5):
        parser.error("--max-order must be between 0 and 5")
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    rows = catalog(args.max_order)
    results = []
    for row in rows:
        print(f"{row['descriptor']:<22} expected={row['expected_dofs']:<4}", end=" ", flush=True)
        result = run_one(row, args.timeout)
        results.append(result)
        observed = result.get("nbdof", "-")
        print(f"observed={str(observed):<4} {result['status']}", flush=True)
    out = {
        "audit_type": "constructor_and_local_dimension_only",
        "source_version_expected": "GetFEM v5.5 tag c242f7fe0989e9a90846c63eaccc6b495d774724",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python_executable": sys.executable,
        "max_order": args.max_order,
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n")
    bad = sum(row["status"] != "constructed_matching_dimension" for row in results)
    print(f"\nConstructor audit: {len(results) - bad}/{len(results)} constructed "
          f"with expected dimensions; {bad} errors/mismatches/timeouts.")
    print(f"JSON report: {args.output}")
    print("WARNING: This is not an H(div)/H(curl), Piola, or MPI validation.")
    return int(args.strict and bad > 0)


if __name__ == "__main__":
    sys.exit(main())
