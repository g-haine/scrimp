# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Original Heat-Wave native DAE incompatibility diagnostic, not a TS PASS.

Run unmodified top-level original Heat-Wave with default L2 projection,
distributed COO, fixed prescribed differential values and the genuine
PETSc algebraic initializer. Expected outcome is an explicit
InconsistentDAEInitialData. Capture actual PETSc/IFunction residual
partition from the production initializer without source or model edits.
"""

import argparse
import json
import math
from pathlib import Path

from mpi4py import MPI
from scrimp.dphs import DPHS
from scrimp.utils.dae_initialization import InconsistentDAEInitialData
from examples.heat_wave import heat_wave_eq


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    comm=MPI.COMM_WORLD
    if comm.size!=1:
        raise RuntimeError("Initial FE/DAE incompatibility diagnostic is MPI1 only")
    if not args.out.parent.is_dir() or args.out.exists():
        raise RuntimeError("Unmodified new output path required")
    captured={}
    original=DPHS.init_step

    def capture_rejection(self):
        try:
            return original(self)
        except InconsistentDAEInitialData as exc:
            captured["reason"]=str(exc)
            captured["report"]=self.initialization_report
            raise

    DPHS.init_step=capture_rejection
    try:
        try:
            heat_wave_eq(
                heat_region=1,wave_region=2,final_time=.004,
                time_step=.001,save_step=.001,factor_solver="mumps",
            )
        except InconsistentDAEInitialData:
            pass
        else:
            raise RuntimeError(
                "Expected incompatible ORIGINAL Gaussian and exterior "
                "homogeneous Dirichlet data to be rejected"
            )
    finally:
        DPHS.init_step=original

    report=captured.get("report")
    if not report or report["status"]!="UNRESOLVED":
        raise RuntimeError("Native DAE rejection missing structured report")
    by_variable=report.get("algebraic_residual_by_variable")
    if not by_variable:
        raise RuntimeError("Native algebraic residual decomposition missing")
    final=float(report["native_final_residual_l2"])
    tol=float(report["feasibility_tolerance"])
    if (final<=tol or not report["differential_initial_exact"]
            or report["pseudo_ts_executed"]):
        raise RuntimeError("Original DAE feasibility semantics violated")
    partitions=sum(v["residual_l2"]**2 for v in by_variable.values())
    if not math.isclose(partitions,final**2,rel_tol=1e-7,abs_tol=1e-18):
        raise RuntimeError("Native block residual does not partition")
    boundary_analytic=5.*math.exp(-4.)
    if boundary_analytic<=0.:
        raise RuntimeError("Analytic Gaussian boundary assertion failed")
    data={
        "schema":"scrimp-v130-heat-wave-original-dae-feasibility-v1",
        "status":"PASS_DETECTED_PHYSICALLY_INCOMPATIBLE_INITIAL_DATA",
        "source_model":"examples/heat_wave.py (unmodified)",
        "numerical_method":"default L2 projection, distributed COO",
        "DAE_strategy":"algebraic_fixed_differential",
        "original_DAE_error":captured["reason"],
        "native_initializer_report":report,
        "analytical_wave_gaussian_at_outer_boundary_x1_y0":
            boundary_analytic,
        "wave_outer_boundary_prescribed_value":0.,
        "full_TS_executed":False,
        "physical_model_qualified":False,
        "source_files_changed":False,
    }
    args.out.write_text(json.dumps(data,indent=2)+"\n",
                        encoding="utf-8")
    print("Native original Heat-Wave DAE incompatibility: "
          "correctly REJECTED",flush=True)
    for row in report["dominant_algebraic_variables"]:
        print(
            f"  algebraic {row['variable']} port={row['port']} "
            f"region={row['region']} n={row['dofs']} "
            f"native_residual={row['residual_l2']:.9g}",
            flush=True,
        )
    print(f"  true ||F_a||={final:.9g}, tolerance={tol:.9g}; "
          f"analytic outer p0(1,0)={boundary_analytic:.9g}",flush=True)
    print("PASS_DETECTED_PHYSICALLY_INCOMPATIBLE_INITIAL_DATA",flush=True)


if __name__=="__main__":
    main()
