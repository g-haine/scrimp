# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             examples/heat.py
- authors:          Giuseppe Ferraro, Ghislain Haine
- date:             22 nov. 2022
- brief:            2D heat equation with Lyapunov Hamiltonian
"""

# Import scrimp
import scrimp as S
import numpy as np
from pathlib import Path

def heat_eq(
    final_time=5.0, time_step=0.01, mesh_size=0.1, *,
    time_stepper="cn", initializer=None, factor_solver=None,
    showcase=False, showcase_output=None, check=False,
    init_profile=None,
):
    """2-D mixed-boundary heat pH-DAE, with the original weak formulation.

    The v1.3 demonstration selects fixed-step Crank--Nicolson with
    strict consistency projection of algebraic unknowns (no pseudo-time
    alteration of the prescribed differential initial state). The
    stricter per-example initialization was necessary for full MPI
    agreement of every coefficient; it does not change SCRIMP's
    global initializer defaults.

    To reproduce the separately qualified historical BDF2 case, call
    heat_eq(time_stepper="bdf", initializer="pseudo", showcase=False).
    These options preserve its original bricks, initial data and
    three Dirichlet + one Neumann/Robin feedback control.

    When showcase=True, the example ITSELF exports the actual GetFEM
    temperature time series and renders scientific PNG/SVG/GIF figures.
    check=True enforces the genuine CN stepwise power balance.
    """
    if time_stepper not in ("cn", "bdf"):
        raise ValueError("time_stepper must be 'cn' or 'bdf'")
    if initializer is None:
        initializer = "algebraic" if time_stepper == "cn" else "pseudo"
    if initializer not in ("algebraic", "pseudo"):
        raise ValueError("initializer must be 'algebraic' or 'pseudo'")
    if init_profile is None:
        # The numerically qualified CN/DAE example must NOT silently
        # select the loose initializer that excites alternating
        # algebraic CN modes with rank-dependent FE coefficients.
        # Preserve the historical BDF2/pseudo path separately.
        init_profile = "strict" if time_stepper == "cn" else "default"
    if init_profile not in ("default", "strict"):
        raise ValueError("init_profile must be 'default' or 'strict'")
    if init_profile == "strict" and initializer != "algebraic":
        raise ValueError("Strict initialization requires algebraic projection")
    if not (np.isfinite(final_time) and np.isfinite(time_step)
            and np.isfinite(mesh_size)
            and 0. < final_time <= 50.
            and 0. < time_step <= final_time and mesh_size > 0.):
        raise ValueError("Invalid heat example mesh/time parameters")
    if (showcase or check) and time_stepper != "cn":
        raise ValueError(
            "Native CN showcase/check requires time_stepper='cn'; "
            "the original BDF2 case remains separately qualified"
        )
    if showcase_output is not None and not showcase:
        raise ValueError("showcase_output requires showcase=True")

    # Init the distributed port-Hamiltonian system
    heat = S.DPHS("real")
    heat.set_output_name("heat",output_dir=showcase_output)

    # Set the domain (using the built-in geometry `Rectangle`)
    # Labels: Omega = 1, Gamma_Bottom = 10, Gamma_Right = 11, Gamma_Top = 12, Gamma_Left = 13
    heat.set_domain(S.Domain("Rectangle", {"L": 2.0, "l": 1.0, "h": float(mesh_size)}))

    # Define the variables and their discretizations and add them to the dphs
    states = [
        S.State("T", "Temperature", "scalar-field"),
    ]
    costates = [
        # Substituted=True indicates that only one variable has to be discretized on this port
        S.CoState("T", "Temperature", states[0], substituted=True)
    ]

    ports = [
        S.Port("Heat flux", "f_Q", "J_Q", "vector-field"),
    ]

    control_ports = [
        S.Control_Port(
            "Boundary control (bottom)",
            "U_B",
            "Temperature",
            "Y_B",
            "- Normal heat flux",
            "scalar-field",
            region=10,
            position="effort",
        ),
        S.Control_Port(
            "Boundary control (right)",
            "U_R",
            "Temperature",
            "Y_R",
            "- Normal heat flux",
            "scalar-field",
            region=11,
            position="effort",
        ),
        S.Control_Port(
            "Boundary control (top)",
            "U_T",
            "Temperature",
            "Y_T",
            "- Normal heat flux",
            "scalar-field",
            region=12,
            position="effort",
        ),
        S.Control_Port(
            "Boundary control (left)",
            "U_L",
            "- Normal heat flux",
            "Y_L",
            "Temperature",
            "scalar-field",
            region=13,
            position="flow",
        ),
    ]

    for state in states:
        heat.add_state(state)
    for costate in costates:
        heat.add_costate(costate)
    for port in ports:
        heat.add_port(port)
    for ctrl_port in control_ports:
        heat.add_control_port(ctrl_port)

    FEMs = [
        S.FEM(states[0].get_name(), 1, FEM="DG"),
        S.FEM(ports[0].get_name(), 2, FEM="CG"),
        S.FEM(control_ports[0].get_name(), 1, FEM="DG"),
        S.FEM(control_ports[1].get_name(), 1, FEM="DG"),
        S.FEM(control_ports[2].get_name(), 1, FEM="DG"),
        S.FEM(control_ports[3].get_name(), 1, FEM="DG"),
    ]
    for FEM in FEMs:
        heat.add_FEM(FEM)

    # Define the physical parameters
    parameters = [
        S.Parameter("rho", "Mass density times heat capacity", "scalar-field", "3.", "T"),
        S.Parameter(
            "Lambda",
            "Heat conductivity",
            "tensor-field",
            "[[1e-2,0.],[0.,1e-2]]",
            "Heat flux",
        ),
    ]
    # Add them to the dphs
    for parameter in parameters:
        heat.add_parameter(parameter)

    # Define the Dirac structure and the constitutive relations block matrices as `Brick`
    bricks = [
        # Add the mass matrices from the left-hand side: the `flow` part of the Dirac structure
        S.Brick("M_T", "T*rho*Test_T", [1], dt=True, position="flow"),
        S.Brick("M_Q", "f_Q.Test_f_Q", [1], position="flow"),
        S.Brick("M_Y_B", "Y_B*Test_Y_B", [10], position="flow"),
        S.Brick("M_Y_R", "Y_R*Test_Y_R", [11], position="flow"),
        S.Brick("M_Y_T", "Y_T*Test_Y_T", [12], position="flow"),
        # Normal trace is imposed by Lagrange multiplier on the left side == the collocated output
        S.Brick("M_Y_L", "U_L*Test_Y_L", [13], position="flow"),
        # Add the matrices from the right-hand side: the `effort` part of the Dirac structure
        S.Brick("D", "-Div(J_Q)*Test_T", [1], position="effort"),
        S.Brick("-D^T", "T*Div(Test_f_Q)", [1], position="effort"),
        S.Brick("B_B", "-U_B*Test_f_Q.Normal", [10], position="effort"),
        S.Brick("B_R", "-U_R*Test_f_Q.Normal", [11], position="effort"),
        S.Brick("B_T", "-U_T*Test_f_Q.Normal", [12], position="effort"),
        # Normal trace is imposed by Lagrange multiplier on the left side == the collocated output
        S.Brick("B_L", "-Y_L*Test_f_Q.Normal", [13], position="effort"),
        S.Brick("C_B", "J_Q.Normal*Test_Y_B", [10], position="effort"),
        S.Brick("C_R", "J_Q.Normal*Test_Y_R", [11], position="effort"),
        S.Brick("C_T", "J_Q.Normal*Test_Y_T", [12], position="effort"),
        S.Brick("C_L", "J_Q.Normal*Test_Y_L", [13], position="effort"),
        ## Define the constitutive relations as getfem `brick`
        # Fourier's law under implicit form - M_e_Q e_Q + CR_Q Q = 0
        S.Brick("-M_J_Q", "-J_Q.Test_J_Q", [1]),
        S.Brick("CR_Q", "f_Q.Lambda.Test_J_Q", [1]),
    ]
    for brick in bricks:
        heat.add_brick(brick)

    # Initialize the problem
    expressions = ["1.", "1.", "1.", "0.2*T"]

    for control_port, expression in zip(control_ports, expressions):
        # Set the control functions (automatic construction of bricks such that -M_u u + f(t) = 0)
        heat.set_control(control_port.get_name(), expression)

    # Set the initial data
    heat.set_initial_value("T", "1. + 2.*np.exp(-50*((x-1)*(x-1)+(y-0.5)*(y-0.5))**2)")

    # The mixed weak forms, co-state substitution and controls above
    # are IDENTICAL in the original BDF2 and the new CN candidate.
    # With constant mass and a quadratic Lyapunov Hamiltonian, CN
    # supports a direct stepwise dissipative port-power identity.
    # Keep the old BDF2/pseudo path EXPLICIT and reproducible.
    scheme = dict(
        init_step_ts_type=initializer,
        t_f=float(final_time),
        ts_type=time_stepper,
        dt=float(time_step),
    )
    if time_stepper == "cn":
        scheme.update(ts_adapt_type="none", dt_save=float(time_step),
                      ts_exact_final_time="matchstep")
    else:
        scheme["ts_bdf_order"] = 2
    if factor_solver is not None:
        scheme["pc_factor_mat_solver_type"] = factor_solver
    if init_profile == "strict":
        # Experimental MPI-parity diagnostic. Pass options HERE,
        # after DPHS construction: DPHS.__init__ clears the global
        # PETSc options database. Setting PETSc.Options before
        # heat_eq() silently loses all such requests.
        scheme.update(
            init_algebraic_max_passes=8,
            init_algebraic_accuracy=1e-13,
            init_algebraic_ksp_rtol=1e-14,
            init_algebraic_ksp_atol=1e-14,
        )
    heat.set_time_scheme(**scheme)

    # Solve
    heat.solve()

    ## Post-processing
    # Set Hamiltonian name
    heat.hamiltonian.set_name("Lyapunov formulation")
    # Define the term
    terms = [
        S.Term("L^2-norm", "0.5*T*rho*T", [1]),
    ]
    # Add them to the Hamiltonian
    for term in terms:
        heat.hamiltonian.add_term(term)

    # This canonical figure and its physical ports remain directly
    # reproducible by running this very example.
    heat.plot_Hamiltonian(save_figure=True, show=False, filename="Hamiltonian_Heat_2D.png")

    if check or showcase:
        observations = measure_heat_cn(heat)
        heat.heat_cn_observations = observations
        if check and not all(observations["gates"].values()):
            raise RuntimeError(
                "Canonical CN heat model failed native GetFEM power gates: "
                + repr(observations["gates"])
            )
        if showcase:
            # All MPI ranks participate in this validated collective.
            heat.export_to_pv("T", t="All")
            from petsc4py import PETSc
            if PETSc.COMM_WORLD.getRank() == 0:
                base = heat._output_root
                heat.heat_visuals = render_heat_visuals(
                    base, observations,
                    heat.get_output_path("mesh") / "T" / "T.pvd",
                )

    return heat


def measure_heat_cn(system):
    """Independently check each native CN step via actual FE/port power.

    Control-port power is U*Y on the true mixed FE ports. Replacing
    this with pointwise T*(J_Q.Normal) is physically WRONG for this
    DG/Dirichlet-multiplier discretization.
    """
    from petsc4py import PETSc

    t = np.asarray(system.solution["t"], dtype=float)
    z = np.asarray(system.solution["z"])
    H = np.asarray(system.get_Hamiltonian(), dtype=float)
    n = len(t)
    if (n < 3 or z.shape[0] != n or H.shape != (n,)
            or not np.all(np.isfinite(z))
            or not np.all(np.isfinite(H))
            or not np.all(np.diff(t) > 0.)
            or system.last_ts_converged_reason <= 0
            or system.last_ts_rejections != 0
            or system.last_ts_snes_failures != 0
            or system.last_ts_steps != n-1):
        raise RuntimeError("Invalid native accepted CN time/FE trajectory")

    if system.time_scheme.getString("ts_type", "").lower() != "cn":
        raise RuntimeError("CN-specific heat balance requested for non-CN TS")
    if system.time_scheme.getString("ts_adapt_type", "").lower() != "none":
        raise RuntimeError("Canonical heat CN must have fixed time steps")
    if not system.powers_computed:
        system.compute_powers()

    # GetFEM evaluates the dissipation on ORIGINAL mid-step FE states.
    d = np.asarray(system.get_quantity(
        "f_Q.Lambda.f_Q", region=1, CN=True
    ), dtype=float)
    port_names = (
        "Boundary control (bottom)", "Boundary control (right)",
        "Boundary control (top)", "Boundary control (left)"
    )
    power = {}
    for name in ("Heat flux", *port_names):
        port = system.ports[name]
        if port.get_power_sampling() != "cn_step":
            raise RuntimeError("CN stepwise power sampling was not used")
        values = np.asarray(port.get_power(), dtype=float)
        if values.shape != t.shape or not np.all(np.isfinite(values)):
            raise RuntimeError("Incomplete true FE port power: " + name)
        power[name] = values
    controls = np.sum(np.stack([power[k] for k in port_names]), axis=0)
    dissipation_port = power["Heat flux"]
    dt = np.diff(t)
    dH = np.diff(H)
    residual = dH + dt*(d[1:]+controls[1:])
    total_scale = max(1., float(np.max(np.abs(H))))
    port_disagreement = float(np.max(np.abs(
        d[1:]-dissipation_port[1:]
    )))
    max_defect = float(np.max(np.abs(residual)))
    cumulative = np.r_[0., np.cumsum(
        dt*(d[1:]+controls[1:])
    )]
    report = {
        "status": "CN_NATIVE_PHYSICS_DIAGNOSTIC",
        "MPI": PETSc.COMM_WORLD.getSize(),
        "steps": int(system.last_ts_steps),
        "DOFs": int(z.shape[1]),
        "H_initial": float(H[0]),
        "H_final": float(H[-1]),
        "max_abs_CN_step_energy_defect": max_defect,
        "max_abs_CN_cumulative_energy_defect": float(
            np.max(np.abs(H-H[0]+cumulative))
        ),
        "max_abs_dissipation_port_disagreement": port_disagreement,
        "gates": {
            "original_native_TS_completed":
                system.last_ts_converged_reason > 0,
            "positive_H": bool(np.min(H) > 0),
            "nonnegative_Fourier_dissipation":
                bool(np.min(d) >= -1e-10),
            "Fourier_constitutive_agreement":
                bool(port_disagreement <= 2e-8*total_scale),
            "CN_step_energy_identity":
                bool(max_defect <= 2e-8*total_scale),
            "CN_cumulative_energy_identity":
                bool(np.max(np.abs(H-H[0]+cumulative))
                     <= 2e-8*total_scale),
        },
    }
    return {"time": t, "H": H, "dissipation": d,
            "control_powers": {key: power[key] for key in port_names},
            "CN_step_defect": np.r_[0., residual],
            "CN_cumulative_defect": H-H[0]+cumulative,
            "cumulative_power": cumulative, "report": report,
            "gates": report["gates"]}


def render_heat_visuals(folder, observed, pvd_file, max_gif_frames=36):
    """Render genuine GetFEM DG1 temperature snapshots and budgets.

    All frames originate from actual saved FE snapshots and native
    GetFEM VTU files. The colormap range is identical across frames,
    and no fictitious spatial or temporal interpolation is used.
    """
    import base64
    import io
    import json
    import struct
    from xml.etree import ElementTree as ET

    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from matplotlib import colors, cm
    from PIL import Image

    folder = Path(folder)
    png=folder/"png"
    svg=folder/"svg"
    gif=folder/"gif"
    metadata=folder/"json"
    for directory in (png,svg,gif,metadata):
        directory.mkdir(parents=True,exist_ok=True)
    t = observed["time"]
    H = observed["H"]
    d = observed["dissipation"]
    controls = observed["control_powers"]
    report = observed["report"]

    def read_array(node, *, field_size=None):
        mapping = {
            "Float32": "<f4", "Float64": "<f8",
            "Int32": "<i4", "Int64": "<i8",
            "UInt8": "u1", "UInt32": "<u4", "UInt64": "<u8",
        }
        if node is None or node.get("type") not in mapping:
            raise ValueError("Invalid original GetFEM VTK DataArray type")
        dtype = np.dtype(mapping[node.get("type")])
        encoding = node.get("format", "ascii")
        if encoding == "ascii":
            arr = np.fromstring(node.text or "", dtype=dtype, sep=" ")
        elif encoding == "binary":
            payload = base64.b64decode(
                "".join((node.text or "").split()), validate=True
            )
            if len(payload) < 4 or (len(payload)-4) % dtype.itemsize:
                raise ValueError("Malformed actual GetFEM VTU binary")
            actual = len(payload)-4
            count = struct.unpack("<I",payload[:4])[0]
            if count != actual:
                # Genuine GetFEM Float32 point fields may encode the
                # number of scalar entries in the 4-byte header.
                if (field_size is None or node.get("type") != "Float32"
                        or struct.unpack("<f",payload[:4])[0]
                           != float(field_size)):
                    raise ValueError("Invalid genuine GetFEM VTU header")
            arr = np.frombuffer(payload,dtype=dtype,offset=4).copy()
        else:
            raise ValueError("Unsupported native GetFEM VTK encoding")
        if not np.all(np.isfinite(arr)):
            raise ValueError("Nonfinite native GetFEM FE visualization")
        return arr

    def read_temperature(filename):
        piece = ET.parse(filename).getroot().find(".//Piece")
        if piece is None:
            raise ValueError("Missing GetFEM native VTK Piece")
        npt = int(piece.get("NumberOfPoints", "0"))
        nc = int(piece.get("NumberOfCells", "0"))
        if npt < 3 or nc < 1:
            raise ValueError("Empty native heat FE triangular mesh")
        coords = read_array(piece.find("./Points/DataArray"))
        if coords.shape != (3*npt,):
            raise ValueError("Malformed genuine GetFEM VTK points")
        pts = coords.reshape(npt,3)
        if np.max(np.abs(pts[:,2]))>1e-7:
            raise ValueError("Heat FE mesh is not planar")
        c = read_array(piece.find("./Cells/DataArray[@Name='connectivity']"))
        o = read_array(piece.find("./Cells/DataArray[@Name='offsets']"))
        k = read_array(piece.find("./Cells/DataArray[@Name='types']"))
        if (c.shape!=(3*nc,) or o.shape!=(nc,) or k.shape!=(nc,)
                or not np.array_equal(o,3*np.arange(1,nc+1))
                or not np.all(k==5)):
            raise ValueError("Original heat DG1 must yield VTK triangles")
        tri = np.asarray(c,dtype=int).reshape(nc,3)
        if np.min(tri)<0 or np.max(tri)>=npt:
            raise ValueError("Invalid heat VTK triangle connectivity")
        field = piece.find("./PointData/DataArray[@Name='T']")
        if field is None:
            raise ValueError("GetFEM heat VTU has no native temperature T")
        values = read_array(field,field_size=npt)
        if values.shape!=(npt,):
            raise ValueError("Unexpected original native DG1 T field")
        return pts[:,:2].astype(float),tri,values.astype(float)

    pvd = ET.parse(pvd_file).getroot().findall(".//DataSet")
    if len(pvd)!=len(t):
        raise ValueError("Canonical heat GetFEM PVD clock/FE length mismatch")
    entries=[]
    for i, frame in enumerate(pvd):
        if abs(float(frame.attrib["timestep"])-float(t[i])) > 1e-9:
            raise ValueError("Heat GetFEM PVD time not on accepted FE clock")
        name=Path(frame.attrib["file"])
        if name.name != str(name) or name.suffix != ".vtu":
            raise ValueError("Unsafe native heat VTK snapshot reference")
        fn=pvd_file.parent/name
        if not fn.is_file() or fn.stat().st_size<100:
            raise ValueError("Missing original heat GetFEM snapshot")
        entries.append(fn)

    selected=np.unique(np.linspace(
        0,len(t)-1,min(len(t),max_gif_frames),dtype=int
    ))
    actual=[]
    geometry=None
    for i in selected:
        pts, tri, temp=read_temperature(entries[int(i)])
        if geometry is None:
            geometry=(pts,tri)
        elif (not np.array_equal(pts,geometry[0])
              or not np.array_equal(tri,geometry[1])):
            raise ValueError("Unexpected heat mesh changes during native run")
        actual.append(temp)
    actual=np.asarray(actual)
    pts,tri=geometry
    vmin=float(np.min(actual))
    vmax=float(np.max(actual))
    if not vmin<vmax:
        raise ValueError("Heat model spatial field is not informative")
    native_cells=pts[tri]
    norm=colors.Normalize(vmin=vmin,vmax=vmax)
    cmap=plt.get_cmap("inferno")

    def field_plot(index,ax):
        # DG1 is discontinuous: triangles carry their own true FE
        # vertex values. Displaying their cellwise mean avoids
        # inventing smooth continuity across DG element faces.
        cells=actual[index][tri].mean(axis=1)
        artist=PolyCollection(
            native_cells, array=cells, cmap=cmap, norm=norm,
            edgecolors="none", linewidths=0
        )
        ax.add_collection(artist)
        ax.set(xlim=(0,2),ylim=(0,1),aspect="equal",
               xlabel="x",ylabel="y",
               title=f"t={t[selected[index]]:.2f}")
        return artist

    # 4 reproducible fixed-color-range actual FE snapshots.
    chosen=sorted(set([0,len(selected)//3,
                       2*len(selected)//3,len(selected)-1]))
    fig,axs=plt.subplots(2,2,figsize=(12.5,6.7),constrained_layout=True)
    for ax,ind in zip(axs.flat,chosen):
        field_plot(ind,ax)
    fig.colorbar(
        cm.ScalarMappable(norm=norm,cmap=cmap),
        ax=axs.ravel().tolist(),shrink=.75,label="T (GetFEM DG1)"
    )
    fig.suptitle("2-D heat diffusion — original GetFEM temperature field")
    fig.savefig(png/"heat_temperature_snapshots.png",dpi=170)
    fig.savefig(svg/"heat_temperature_snapshots.svg")
    plt.close(fig)

    # Physical balance: CN step work, NOT an approximated BDF2 integral.
    cumulative_d=np.r_[0.,np.cumsum(np.diff(t)*d[1:])]
    cumulative_c=np.r_[0.,np.cumsum(np.diff(t)*
        np.sum(np.stack([p[1:] for p in controls.values()]),axis=0))]
    fig,axes=plt.subplots(2,2,figsize=(13.,8.),constrained_layout=True)
    axes[0,0].plot(t,H,label="H(t) · genuine GetFEM")
    axes[0,0].set(ylabel="Hamiltonian H",title="Lyapunov energy")
    axes[0,0].legend()
    axes[0,1].plot(t,cumulative_d,label="Fourier dissipation")
    axes[0,1].plot(t,cumulative_c,label="Sum of boundary control powers")
    axes[0,1].set(ylabel="Cumulative power",title="Physical exchanges")
    axes[0,1].legend()
    axes[1,0].plot(t,H-H[0],label="H(t) − H(0)")
    axes[1,0].plot(t,-cumulative_d-cumulative_c,
                   "--",label="− D − port work (CN)")
    axes[1,0].set(ylabel="ΔH",title="Stepwise physical power balance")
    axes[1,0].legend()
    axes[1,1].plot(t,observed["CN_cumulative_defect"],
                   label="H(t)−H(0)+∫(D+P) CN steps")
    axes[1,1].set(ylabel="Balance defect",title="Discrete CN closure")
    axes[1,1].legend()
    for ax in axes.flat:
        ax.set_xlabel("Time")
        ax.grid(True,alpha=.25)
    fig.suptitle("2-D heat · real GetFEM dissipation and control ports")
    fig.savefig(png/"heat_energy_and_ports.png",dpi=170)
    fig.savefig(svg/"heat_energy_and_ports.svg")
    plt.close(fig)

    # One native DG1 field per frame, with common color limits and no
    # temporal interpolation. The original complete VTK/PVD stays
    # available in outputs/heat/mesh/T for ParaView.
    gif_frames=[]
    for k in range(len(selected)):
        fig,ax=plt.subplots(figsize=(7.2,3.9))
        field_plot(k,ax)
        fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax,
                     shrink=.8,label="T (native DG1)")
        fig.suptitle("SCRIMP · heat diffusion · actual FE snapshot")
        fig.tight_layout()
        buffer=io.BytesIO()
        fig.savefig(buffer,format="png",dpi=90)
        plt.close(fig)
        buffer.seek(0)
        with Image.open(buffer) as image:
            gif_frames.append(image.convert("RGB").copy())
    gif_frames[0].save(
        gif/"heat_temperature_evolution.gif",
        save_all=True,append_images=gif_frames[1:],
        duration=170,loop=0,optimize=False,
    )

    metadata={
        **report,
        "visual_status":"GENERATED_FROM_ACTUAL_GETFEM_NATIVE_FE",
        "temperature_png":"heat_temperature_snapshots.png",
        "energy_png":"heat_energy_and_ports.png",
        "temperature_gif":"heat_temperature_evolution.gif",
        "genuine_VTK_PVD":str(pvd_file),
        "VTK_frames":len(entries),
        "GIF_original_FE_frames":len(selected),
        "same_temperature_colormap_all_frames":[vmin,vmax],
        "native_DG1_cellwise_graphics_only":True,
        "no_spatial_or_time_interpolation":True,
    }
    (metadata/"manifest.json").write_text(
        json.dumps(metadata,indent=2)+"\n",encoding="utf-8"
    )
    return metadata


def main():
    import argparse
    import json
    from petsc4py import PETSc

    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--time-stepper",choices=("cn","bdf"),default="cn")
    parser.add_argument("--initializer",choices=("algebraic","pseudo"),default=None)
    parser.add_argument("--init-profile",choices=("default","strict"),
                        default=None,
                        help="Algebraic precision profile: strict for the "
                             "CN demo by default, default for legacy BDF2")
    parser.add_argument("--solver",choices=("mumps","superlu","superlu_dist"),
                        default="mumps")
    parser.add_argument("--mesh-size",type=float,default=.1)
    parser.add_argument("--final-time",type=float,default=5.)
    parser.add_argument("--time-step",type=float,default=.01)
    parser.add_argument("--showcase",dest="showcase",action="store_true")
    parser.add_argument("--no-showcase",dest="showcase",action="store_false")
    parser.set_defaults(showcase=True)
    parser.add_argument("--showcase-output",type=Path,default=None)
    parser.add_argument("--check",action="store_true",
                        help="Enforce real CN stepwise native FE power identity")
    parser.add_argument("--output",type=Path,default=None,
                        help="Optional output directory for independent JSON/NPZ traces")
    args=parser.parse_args()
    effective_init_profile = args.init_profile or (
        "strict" if args.time_stepper == "cn" else "default"
    )
    if args.solver == "superlu" and PETSc.COMM_WORLD.getSize()>1:
        parser.error("Sequential SuperLU is incompatible with MPI>1")
    if (args.showcase or args.check) and args.time_stepper != "cn":
        parser.error("Legacy BDF2 runs require --no-showcase and no --check")
    if args.showcase_output is not None and not args.showcase:
        parser.error("--showcase-output requires showcase")

    system=heat_eq(
        final_time=args.final_time,time_step=args.time_step,
        mesh_size=args.mesh_size,time_stepper=args.time_stepper,
        initializer=args.initializer,factor_solver=args.solver,
        showcase=args.showcase,showcase_output=args.showcase_output,
        check=args.check,init_profile=effective_init_profile,
    )
    report=getattr(system,"heat_cn_observations",None)
    if args.output is not None:
        args.output.mkdir(parents=True,exist_ok=True)
        rank=PETSc.COMM_WORLD.getRank()
        import hashlib
        z=np.asarray(system.solution["z"])
        data={
            "time":np.asarray(system.solution["t"],dtype=float),
            "states":z,
            "H":np.asarray(system.get_Hamiltonian(),dtype=float),
        }
        if report is not None:
            data.update({
                "Fourier_dissipation":report["dissipation"],
                "CN_step_defect":report["CN_step_defect"],
                "CN_cumulative_defect":report["CN_cumulative_defect"],
                **{"power_"+str(k):v for k,v in
                   enumerate(report["control_powers"].values())}
            })
        np.savez_compressed(args.output/f"rank{rank}.npz",**data)
        metadata={
            "status":("PASS_CN_HEAT_NATIVE_PHYSICS"
                      if report and all(report["gates"].values())
                      else "BDF2_LEGACY_EXECUTION_ONLY"),
            "MPI":PETSc.COMM_WORLD.getSize(),
            "rank":rank,
            "time_stepper":args.time_stepper,
            "initializer":(args.initializer or
                           ("algebraic" if args.time_stepper=="cn" else "pseudo")),
            "init_profile":effective_init_profile,
            "mesh_size":args.mesh_size,
            "time_step":args.time_step,
            "final_time":args.final_time,
            "steps":int(system.last_ts_steps),
            "DOFs":int(z.shape[1]),
            "initialization_report":getattr(system,"initialization_report",None),
            "init_profile_applied":({
                "init_algebraic_max_passes":
                    system.time_scheme.getString("init_algebraic_max_passes","3"),
                "init_algebraic_accuracy":
                    system.time_scheme.getString("init_algebraic_accuracy","1e-11"),
                "init_algebraic_ksp_rtol":
                    system.time_scheme.getString("init_algebraic_ksp_rtol","1e-13"),
                "init_algebraic_ksp_atol":
                    system.time_scheme.getString("init_algebraic_ksp_atol","1e-13"),
            } if effective_init_profile == "strict" else None),
            "state_SHA256":hashlib.sha256(
                np.ascontiguousarray(z).tobytes()
            ).hexdigest(),
            "physics":report["report"] if report else None,
            "visuals":getattr(system,"heat_visuals",None),
        }
        (args.output/f"rank{rank}.json").write_text(
            json.dumps(metadata,indent=2,default=str)+"\n",encoding="utf-8"
        )
    if PETSc.COMM_WORLD.getRank()==0 and report is not None:
        print(json.dumps(report["report"],indent=2),flush=True)
    return system


if __name__ == "__main__":
    main()
