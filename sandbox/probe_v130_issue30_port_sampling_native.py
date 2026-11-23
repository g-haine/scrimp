# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""#30 genuine GetFEM/PETSc native Port.compute + DPHS.plot_powers smoke.

Two original linear wave pH-DAE runs: CN and BDF2, same PDE and
time-dependent boundary law. BDF nodal power is ONLY an approximate
postprocessing observable, never called exact BDF discrete work.
No model, solver or constitutive weak-form modification.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from probe_v130_checkpoint_linear_dae import (
    DT, build_linear_wave, original_IFunction, get_indices,
)
from scrimp.utils.power_integration import integrate_port_power


FINAL=.024
SCHEMES=("cn","bdf")


def run():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--scheme",choices=SCHEMES,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    comm=PETSc.COMM_WORLD.tompi4py()
    system=build_linear_wave(FINAL)
    if args.scheme=="bdf":
        # Only changes the PETSc integrator. The original FEM,
        # time-dependent boundary law and exact rest state stay fixed.
        system.set_time_scheme(
            ts_type="bdf",ts_bdf_order=2,ts_adapt_type="none",
            t_0=0.,t_f=FINAL,dt=DT,dt_save=DT*(1.-1e-9),
            init_step=False,
        )
    if system.time_scheme.getString("ts_type").lower()!=args.scheme:
        raise RuntimeError("Requested scheme not installed in production DPHS")
    system.solve()
    times=np.asarray(system.solution["t"],dtype=float)
    states=np.asarray(system.solution["z"],dtype=float)
    steps=int(round(FINAL/DT))
    if (system.last_ts_converged_reason<=0
            or system.last_ts_rejections
            or system.last_ts_snes_failures
            or system.last_ts_steps!=steps
            or abs(float(system.last_ts_time)-FINAL)>1e-12
            or times.shape!=(steps+1,)
            or states.ndim!=2 or states.shape[0]!=len(times)
            or not np.allclose(times,np.linspace(0.,FINAL,steps+1),
                               rtol=0.,atol=1e-12)
            or not np.all(np.isfinite(states))):
        raise RuntimeError("Original physical GetFEM/PETSc wave TS failed")

    dynamic,algebraic=get_indices(system)
    zeros=np.zeros(states.shape[1])
    max_fa=0.
    for k,t in enumerate(times):
        residual=original_IFunction(system,comm,float(t),states[k],zeros)
        max_fa=max(max_fa,float(np.linalg.norm(residual[algebraic])))
    if max_fa>=1e-8:
        raise RuntimeError(f"Original physical DAE algebraic residual {max_fa}")

    # Forces the ORIGINAL GetFEM assembly and public SCRIMP method.
    H=np.asarray(system.get_Hamiltonian(),dtype=float)
    system.compute_powers()
    sampling="cn_step" if args.scheme=="cn" else "instantaneous"
    convention=("cn_step" if args.scheme=="cn"
                else "instantaneous_trapezoid")
    ports={}
    integrations={}
    for port in system.ports.values():
        if port.get_algebraic() and not port.get_substituted():
            name=port.get_name()
            if port.get_power_sampling()!=sampling:
                raise RuntimeError(
                    f"Actual original GetFEM Port.compute sampling wrong: {name}"
                )
            powers=np.asarray(port.get_power(),dtype=float)
            if powers.shape!=times.shape or not np.all(np.isfinite(powers)):
                raise RuntimeError("Original public port-power samples malformed")
            ports[name]=powers
            integrations[name]=integrate_port_power(
                times,powers,convention=convention
            )
    if not ports or not np.all(np.isfinite(H)):
        raise RuntimeError("No original physical wave port powers or Hamiltonian")
    amplitude=max(float(np.max(abs(x))) for x in ports.values())
    if amplitude<=1e-9:
        raise RuntimeError("Original excited physical wave has no port activity")

    # Each MPI rank must see the SAME generic power for the same
    # GetFEM model. Never infer rank parity from the rank0 PNG.
    rank_data={name:values.tolist() for name,values in ports.items()}
    copies=comm.allgather(rank_data)
    if any(set(d)!=set(ports) for d in copies):
        raise RuntimeError("MPI port sets differ across actual ranks")
    max_rank_diff=0.
    for d in copies:
        for name,values in ports.items():
            max_rank_diff=max(
                max_rank_diff,
                float(np.max(abs(np.asarray(d[name])-values)))
            )
    if max_rank_diff>=1e-8:
        raise RuntimeError(f"Rank-local GetFEM powers differ: {max_rank_diff}")

    net=np.sum(np.stack(list(ports.values())),axis=0)
    cn_max_defect=None
    if args.scheme=="cn":
        per_step=np.diff(H)+np.diff(times)*net[1:]
        cn_max_defect=float(np.max(abs(per_step)))
        if cn_max_defect>=1e-11:
            raise RuntimeError("Real original wave CN power identity failed")
    # For BDF, there is NO claim that integrating nodal values by
    # trapezoids equals the actual implicit multistep work.

    if comm.rank==0:
        if args.out.exists() and any(args.out.iterdir()):
            raise FileExistsError("Refusing to replace public plotting evidence")
        args.out.mkdir(parents=True,exist_ok=True)
        fig,ax=plt.subplots(figsize=(9,5))
    else:
        fig=ax=None

    # This is the exact PUBLIC plotting method under qualification.
    system.plot_powers(ax=ax,HamTot=H)

    if comm.rank==0:
        by_label={line.get_label():line for line in ax.get_lines()}
        nports=0
        expected_balance=H.copy()
        for port in system.ports.values():
            if port.get_algebraic() and not port.get_substituted():
                name=port.get_name()
                sign=-1. if port.get_dissipative() else 1.
                label=(
                    name+" ("+(
                        "CN step integral" if sampling=="cn_step" else
                        "nodal trapezoid (approximate)"
                    )+")"
                )
                if label not in by_label:
                    raise RuntimeError("Public Matplotlib missing labeled port "+name)
                actual=np.asarray(by_label[label].get_ydata(),dtype=float)
                if not np.allclose(actual,sign*integrations[name],
                                   rtol=0.,atol=1e-12):
                    raise RuntimeError("Public plotting port-integral mismatch "+name)
                expected_balance-=sign*integrations[name]
                nports+=1
        candidates=[v for k,v in by_label.items()
                    if k.startswith("Balance (legacy port signs")]
        if len(candidates)!=1 or nports!=len(ports):
            raise RuntimeError("Public Balance curve or number of ports wrong")
        if not np.allclose(candidates[0].get_ydata(),expected_balance,
                           rtol=0.,atol=1e-12):
            raise RuntimeError("Public model-dependent balance curve mismatch")
        ax.set_xlabel("Original physical time")
        fig.tight_layout()
        fig.savefig(args.out/"public_plot_powers.png",dpi=150)
        plt.close(fig)
        if (args.out/"public_plot_powers.png").stat().st_size<1000:
            raise RuntimeError("Public matplotlib plot was empty")
        np.savez_compressed(
            args.out/"original_physical_wave_power.npz",
            times=times,states=states,H=H,
            port_names=np.asarray(list(ports)),
            raw_ports=np.stack(list(ports.values())),
            signed_port_integrals=np.stack(list(integrations.values())),
        )
        report={
            "schema":"scrimp-v130-issue30-port-generic-native-smoke-v1",
            "status":"NATIVE_PUBLIC_PORT_COMPUTE_AND_PLOT_POWERS_PASS",
            "scheme":args.scheme,
            "sampling":sampling,"integration_convention":convention,
            "mpi_ranks":comm.size,"physical_model":"original_wave_2D_linear_DAE",
            "steps":steps,"time_final":float(times[-1]),
            "finite_element_dofs":int(states.shape[1]),
            "getfem_algebraic_residual_max_l2":max_fa,
            "original_ports_checked":list(ports),
            "maximum_original_generic_port_activity":amplitude,
            "maximum_rank_local_port_difference":max_rank_diff,
            "maximum_CN_per_step_signed_energy_defect":cn_max_defect,
            "BDF_nodal_plot_is_approximate_not_discrete_quadrature":
                args.scheme=="bdf",
            "physics_or_solver_modified":False,
        }
        (args.out/"native_public_power_report.json").write_text(
            json.dumps(report,indent=2)+"\n",encoding="utf-8"
        )
        print(
            "NATIVE PUBLIC PLOT POWER PASS "
            f"scheme={args.scheme} MPI={comm.size} DOFs={states.shape[1]} "
            f"sampling={sampling} Fa={max_fa:.9g} "
            f"rank_port_diff={max_rank_diff:.9g} "
            f"CN_step_defect={cn_max_defect}",
            flush=True,
        )
    comm.Barrier()


if __name__=="__main__":
    run()
