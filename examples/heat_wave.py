# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             sandbox/heat_hw.py
- authors:          Giuseppe Ferraro, Ghislain Haine
- date:             15 dec. 2022
- brief:            a 2D coupled heat-wave system
"""

# Import scrimp
import scrimp as S
import numpy as np

from scrimp.utils.fe_initialization import enforce_zero_dirichlet_initial_values


def _read_native_vtk_data_array(element):
    """Read a real GetFEM ASCII or inline-binary VTU numeric DataArray.

    This viewer deliberately uses the *exported GetFEM field* and native
    triangles, never a synthetic wave or a guessed triangulation.
    """
    import base64
    import struct

    if element is None:
        raise ValueError("Missing native GetFEM VTK DataArray")
    name = element.get("type", "")
    types = {
        "Float32": "<f4", "Float64": "<f8",
        "Int32": "<i4", "Int64": "<i8",
        "UInt8": "u1", "UInt32": "<u4", "UInt64": "<u8",
    }
    if name not in types:
        raise ValueError("Unsupported native VTK array type: " + name)
    dtype = np.dtype(types[name])
    form = element.get("format", "ascii")
    if form == "ascii":
        data = np.fromstring(element.text or "", sep=" ", dtype=dtype)
    elif form == "binary":
        raw = base64.b64decode((element.text or "").strip(), validate=True)
        if len(raw) < 4:
            raise ValueError("Truncated native binary VTU")
        length = struct.unpack("<I", raw[:4])[0]
        # GetFEM's PointData is sometimes prefixed with a component
        # counter rather than a byte count. Validate dimensions at
        # the VTK Piece layer independently.
        if length != len(raw)-4 and (len(raw)-4) % dtype.itemsize:
            raise ValueError("Malformed binary VTU payload")
        data = np.frombuffer(raw[4:], dtype=dtype).copy()
    else:
        raise ValueError("Unsupported GetFEM VTU data format: " + form)
    if not np.all(np.isfinite(data)):
        raise ValueError("Nonfinite native GetFEM visualization data")
    return data


def _native_wave_vtk_frame(pvd_dir, frame_file):
    """Read the actual scalar p field on genuine 2-D VTK triangles."""
    from xml.etree import ElementTree as ET

    filename = Path(frame_file)
    if (filename.name != str(filename) or filename.suffix != ".vtu"):
        raise ValueError("Invalid native VTK snapshot path")
    root = ET.parse(pvd_dir / filename).getroot()
    piece = root.find(".//Piece")
    if piece is None:
        raise ValueError("Missing GetFEM 2-D VTK Piece")
    npts = int(piece.get("NumberOfPoints", "0"))
    ncells = int(piece.get("NumberOfCells", "0"))
    if npts <= 0 or ncells <= 0:
        raise ValueError("Empty 2-D GetFEM VTK mesh")
    points = _read_native_vtk_data_array(
        piece.find("./Points/DataArray")
    )
    if points.size != 3*npts:
        raise ValueError("Wrong GetFEM VTU coordinate dimensions")
    points = points.reshape(npts, 3)
    conn = _read_native_vtk_data_array(
        piece.find("./Cells/DataArray[@Name='connectivity']")
    )
    offsets = _read_native_vtk_data_array(
        piece.find("./Cells/DataArray[@Name='offsets']")
    )
    types = _read_native_vtk_data_array(
        piece.find("./Cells/DataArray[@Name='types']")
    )
    if not (conn.size==3*ncells and offsets.size==ncells
            and np.array_equal(offsets,3*np.arange(1,ncells+1))
            and types.size==ncells and np.all(types==5)):
        raise ValueError("Expected authentic triangular 2-D GetFEM mesh")
    triangle = np.asarray(conn.reshape(ncells,3),dtype=int)
    if np.any(triangle<0) or np.any(triangle>=npts):
        raise ValueError("Invalid GetFEM triangle connectivity")
    field = piece.find("./PointData/DataArray[@Name='p']")
    if field is None or int(field.get("NumberOfComponents","1"))!=1:
        raise ValueError("Missing native scalar velocity p")
    values = _read_native_vtk_data_array(field)
    if values.shape!=(npts,):
        raise ValueError("Wrong native scalar wave field length")
    if np.max(np.abs(points[:,2])) > 1e-10:
        raise ValueError("Expected planar GetFEM 2-D geometry")
    return points[:,:2],triangle,values


def render_heat_wave_long_showcase(times, H, *, profile, wave_region):
    """Create genuine log-log and annular-wave visuals from THIS run.

    Only wave_region=2 can provide the annular gallery. VTK nodal
    velocity plots are qualitative: they are not a proof of a
    whispering-gallery eigenmode nor a numerical energy conservation
    gate. No time interpolation, manufactured data or new physics.
    """
    from pathlib import Path
    from xml.etree import ElementTree as ET
    from matplotlib import pyplot as plt
    from matplotlib import tri as mtri
    from PIL import Image
    import io
    import json

    if wave_region != 2:
        raise ValueError("Gallery visualization requires an exterior wave annulus")
    times = np.asarray(times,dtype=float).reshape(-1)
    H = np.asarray(H,dtype=float).reshape(-1)
    if (len(times)<5 or H.shape!=times.shape
            or not np.all(np.isfinite(times))
            or not np.all(np.isfinite(H))
            or not np.all(np.diff(times)>0)
            or H[0]<=0):
        raise ValueError("Incorrect true long-time physical Hamiltonian")
    base = Path.cwd()/"outputs"
    gallery = base/"showcase"/"heat_wave_annulus"
    gallery.mkdir(parents=True,exist_ok=True)

    catalog = base/"pv"/"p"/"p.pvd"
    pvd = ET.parse(catalog).getroot().findall(".//DataSet")
    if len(pvd)!=len(times):
        raise ValueError("GetFEM p VTU frames do not match saved FE times")
    all_times = np.asarray(
        [float(frame.attrib["timestep"]) for frame in pvd],
        dtype=float,
    )
    if not np.allclose(all_times,times,atol=1e-9,rtol=1e-9):
        raise ValueError("GetFEM VTK clock differs from physical FE trajectory")

    # Select actual saved CN instants without interpolating. A long
    # history may have hundreds of snapshots, but a readable GIF needs
    # only 48 representative GetFEM frames; ParaView retains ALL frames.
    selected = np.unique(np.linspace(
        0,len(times)-1,min(len(times),48),dtype=int,
    ))
    positions = None
    triangles = None
    velocities = []
    for i in selected:
        xy,cells,p = _native_wave_vtk_frame(
            catalog.parent,pvd[int(i)].attrib["file"],
        )
        if positions is None:
            positions,triangles = xy,cells
        elif (not np.array_equal(xy,positions)
                or not np.array_equal(cells,triangles)):
            raise ValueError("Wave VTK geometry changed within same trajectory")
        velocities.append(p)
    velocities = np.asarray(velocities)
    if np.max(np.abs(velocities))<=0:
        raise RuntimeError("No nontrivial native annular wave velocity")
    centroid = positions[triangles].mean(axis=1)
    r = np.linalg.norm(centroid,axis=1)
    # Exact Concentric annulus represented by original 2-D GetFEM cells;
    # the centroid cut is only a visualization filter.
    mask = (r>=0.6-1e-5) & (r<=1.+1e-5)
    if np.count_nonzero(mask)<10:
        raise ValueError("No annular GetFEM triangles in native VTK mesh")
    plot_triangles=triangles[mask]
    tri=mtri.Triangulation(
        positions[:,0],positions[:,1],plot_triangles,
    )
    vmax=float(np.max(np.abs(velocities)))
    scale=max(vmax,1e-14)
    theta=np.linspace(0.,2.*np.pi,300)
    def plot_wave(index):
        fig,ax=plt.subplots(figsize=(6.8,6.4))
        painting=ax.tripcolor(
            tri,velocities[index],shading="gouraud",
            cmap="RdBu_r",vmin=-scale,vmax=scale,
        )
        for radius in (0.6,1.):
            ax.plot(radius*np.cos(theta),radius*np.sin(theta),
                    color="black",linewidth=0.65,alpha=0.7)
        fig.colorbar(painting,ax=ax,shrink=0.8,label="p (native GetFEM)")
        ax.set(aspect="equal",xlabel="x",ylabel="y",xlim=(-1.05,1.05),
               ylim=(-1.05,1.05))
        ax.set_title(
            f"Heat–Wave: velocity in outer annulus · t={times[selected[index]]:.3f}\n"
            f"{profile} · native FE/VTK (not an eigenmode proof)",
        )
        fig.tight_layout()
        return fig

    # Log-log uses ONLY actual FE Hamiltonian values and t>0.
    positive=(times>0.) & (H>0.)
    if not np.count_nonzero(positive)>=3:
        raise ValueError("Not enough positive real FE Hamiltonian samples")
    fig,ax=plt.subplots(figsize=(8.6,5.2))
    ax.loglog(times[positive],H[positive],linewidth=1.8,
              label="H(t), native GetFEM quadratures")
    ax.set(xlabel="Time",ylabel="Total Hamiltonian H(t)",
           title=f"Long-time Heat–Wave Hamiltonian · {profile}")
    ax.grid(True,which="both",alpha=.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(gallery/"heat_wave_H_long_loglog.png",dpi=170)
    fig.savefig(gallery/"heat_wave_H_long_loglog.svg")
    plt.close(fig)

    snapshots = sorted(set((0,len(selected)//3,
                            2*len(selected)//3,len(selected)-1)))
    for index in snapshots:
        fig=plot_wave(index)
        fig.savefig(
            gallery/f"wave_annulus_t{times[selected[index]]:.3f}.png",
            dpi=150,
        )
        plt.close(fig)

    frames=[]
    for i in range(len(selected)):
        fig=plot_wave(i)
        buffer=io.BytesIO()
        fig.savefig(buffer,format="png",dpi=96)
        plt.close(fig)
        buffer.seek(0)
        with Image.open(buffer) as im:
            frames.append(im.convert("RGB").copy())
    frames[0].save(
        gallery/"wave_annulus_long.gif",
        save_all=True,append_images=frames[1:],
        duration=180,loop=0,optimize=False,
    )
    # Kinetic concentration diagnostic only: CG1 p^2 exactly integrated
    # per triangle, but the shell r>=0.82 is approximated by cell
    # centroid inclusion. It is NOT a full (p^2+|q|^2) FE energy test.
    area = 0.5*np.abs(
        (positions[triangles[:,1],0]-positions[triangles[:,0],0])*
        (positions[triangles[:,2],1]-positions[triangles[:,0],1])-
        (positions[triangles[:,1],1]-positions[triangles[:,0],1])*
        (positions[triangles[:,2],0]-positions[triangles[:,0],0])
    )
    shell = mask & (r>=0.82)
    fractions=[]
    for field in velocities:
        x=field[triangles]
        kinetic=area/6.*(
            (x*x).sum(axis=1)+x[:,0]*x[:,1]+
            x[:,0]*x[:,2]+x[:,1]*x[:,2]
        )
        denom=float(np.sum(kinetic[mask]))
        fractions.append(
            None if denom<=1e-22
            else float(np.sum(kinetic[shell])/denom)
        )
    valid=np.array([v is not None for v in fractions],dtype=bool)
    if np.any(valid):
        fig,ax=plt.subplots(figsize=(8.5,4.6))
        ax.plot(times[selected][valid],
                np.asarray(fractions,dtype=object)[valid].astype(float),
                marker=".",linewidth=1.1)
        ax.set(xlabel="Time",ylabel="Fraction of CG1 kinetic proxy in r>=0.82",
               title="Outer-shell wave localization (diagnostic only)",
               ylim=(0,1.05))
        ax.grid(alpha=.25)
        fig.tight_layout()
        fig.savefig(gallery/"wave_outer_shell_kinetic_proxy.png",dpi=155)
        plt.close(fig)
    # Reproducible native energy series, not coordinates digitized from
    # a graphic. These are the actual Hamiltonian quadratures at accepted
    # saved FE states.
    np.savez_compressed(
        gallery/"heat_wave_FE_H_long.npz",
        t=np.ascontiguousarray(times,dtype=np.float64),
        H=np.ascontiguousarray(H,dtype=np.float64),
    )
    result={
        "schema":"scrimp-heat-wave-long-native-showcase-v1",
        "status":"VISUALS_GENERATED_FROM_NATIVE_GETFEM_NOT_WG_MODE_CERTIFIED",
        "initial_profile":profile,
        "source":"examples/heat_wave.py same-run GetFEM p.pvd/VTU + FE Hamiltonian",
        "timesteps_saved":int(len(times)),
        "gif_frames_real":int(len(selected)),
        "points":int(len(positions)),
        "annulus_triangles":int(np.count_nonzero(mask)),
        "outer_shell_threshold_radius":0.82,
        "kinetic_proxy_shell_fraction":fractions,
        "selected_true_FE_times":times[selected].tolist(),
        "initial_H":float(H[0]),"final_H":float(H[-1]),
        "no_time_interpolation":True,
        "full_field_energy_localization_certified":False,
        "whispering_gallery_eigenmode_certified":False,
        "note":"p nodal CG1 triangle values only; shell uses centroid rule. "
               "The animation is a visualization, not a modal/eigenvalue proof.",
    }
    (gallery/"provenance.json").write_text(
        json.dumps(result,indent=2)+"\n",encoding="utf-8",
    )
    return result


def heat_wave_eq(
    heat_region=1, wave_region=2, *, factor_solver=None,
    final_time=2.0, time_step=0.001, save_step=0.05,
    initial_wave_boundary="clamp",
    initial_wave_profile="gaussian", long_showcase=False,
):
    """A physically oriented, energy-dissipative coupled heat-wave example.

    Reference: Haine, Matignon & Monteghetti, Mathematical Reports (2022),
    eqs. (2.1)-(2.2), (3.15), (4.3)-(4.5). This example is NOT the paper's
    numerical experiment: it keeps the original SCRIMP heat div-div and
    wave grad-grad weak forms, choosing compatible mixed FEM spaces.

    The original built-in Concentric geometry has disk=1, annulus=2,
    interface=10, outer boundary=20. CN and the 0-input external
    Dirichlet condition are used to expose the heat dissipation and
    closed interconnection in an interpretable power balance.

    Interface normals always mean each subdomain's outward normal.
    The conservative two-sided conforming interface compiler is
    enabled by default in the v1.3.0 release-preparation source.
    Physical labels, FEM declarations and weak forms are unchanged.
    Ambiguous or unsupported interfaces fail explicitly; use
    SCRIMP_AUTO_INTERFACE=0 only to recover legacy manual assembly
    if your model was already oriented by hand.

    Args:
        heat_region: domain label of the heat equation (1 or 2).
        wave_region: domain label of the wave equation (the other one).
        factor_solver: optional PETSc LU solver (default: SCRIMP MUMPS).
        final_time: physical end time.
        time_step: fixed Crank–Nicolson timestep.
        save_step: accepted trajectory sampling interval.
        initial_wave_profile: "gaussian" retains the historical heat/wave
            condition. "gallery_packet" is a smooth high-angular-wavenumber
            wave excitation concentrated near the outer radius, with
            zero initial heat and stress. It is NOT a certified cavity
            eigenmode and is available only for a wave annulus.
        long_showcase: produce an opt-in real-GetFEM annular-wave PNG/GIF
            showcase after successful physical integration and VTK export.
        initial_wave_boundary: compatible wave velocity initialization
            when the exterior boundary belongs to the wave region.
            "clamp" (default) sets exactly the exterior CG1 degrees
            of freedom to zero, leaving interior L2-projected
            coefficients unchanged. "constrained_l2" computes the
            closest zero-trace field in the physical GetFEM L2 norm
            (MPI1 only). "unconstrained" reproduces the historical,
            generally DAE-incompatible Gaussian for diagnostics only.
            Neither valid method relaxes algebraic DAE constraints.
    """
    
    # Init the distributed port-Hamiltonian system
    hw = S.DPHS("real")

    # Set the domain (using the built-in geometry `Concentric`)
    # Labels: Disk = 1, Annulus = 2, Interface = 10, Boundary = 20
    omega = S.Domain("Concentric", {"R": 1.0, "r": 0.6, "h": 0.05})
    
    # And add it to the dphs
    hw.set_domain(omega)
    
    # Define the states and costates, needs the heat and wave region's labels
    states = [
        S.State("T", "Temperature", "scalar-field", region=heat_region),
        S.State("p", "Velocity", "scalar-field", region=wave_region),
        S.State("q", "Stress", "vector-field", region=wave_region),
    ]
    # Use of the `substituted=True` keyword to get the co-energy formulation
    costates = [
        S.CoState("T", "Temperature", states[0], substituted=True),
        S.CoState("p", "Velocity", states[1], substituted=True),
        S.CoState("q", "Stress", states[2], substituted=True),
    ]
    
    # Add them to the dphs
    for state in states:
        hw.add_state(state)
    for costate in costates:
        hw.add_costate(costate)

    # Define the algebraic port
    ports = [
        S.Port("Heat flux", "e_Q", "e_Q", "vector-field", substituted=True, region=heat_region),
    ]
    
    # Add it to the dphs
    for port in ports:
        hw.add_port(port)

    if {heat_region, wave_region} != {1, 2}:
        raise ValueError("Concentric requires heat/wave regions 1 and 2")

    # The outer boundary is either wave velocity (Dirichlet) or heat
    # temperature (Dirichlet). Its port causality depends on the
    # adjacent subsystem, NOT on the user-selected LU factorization.
    outer_is_wave = (wave_region == 2)
    # On the interface, Normal is the outward normal of the bulk
    # subdomain owning the tested variable. The automatic conforming
    # interface compiler selects that side, independently of whether
    # heat occupies the disk or the annulus. The weak transmission law
    # below therefore needs NO correction based on region numbering.

    # Define the control ports
    control_ports = [
        S.Control_Port(
            "Interface Heat", 
            "U_T", 
            "Heat flux", 
            "Y_T", 
            "Temperature",
            "scalar-field",
            region=10, 
            position="effort"
        ),
        S.Control_Port(
            "Interface Wave", 
            "U_w", 
            "Velocity", 
            "Y_w", 
            "Velocity", 
            "scalar-field", 
            region=10, 
            position="effort"
        ),
        # Homogeneous Dirichlet: wave velocity=0 (wave annulus)
        # or temperature=0 (heat annulus). Normal traction/flux is
        # the corresponding Lagrange-multiplier output.
        S.Control_Port(
            "Boundary",
            "U_bnd",
            "Velocity" if outer_is_wave else "Temperature",
            "Y_bnd",
            "Normal force" if outer_is_wave else "Normal heat flux",
            "scalar-field",
            region=20,
            position="flow" if outer_is_wave else "effort",
        ),
    ]

    # Add them to the dphs
    for ctrl_port in control_ports:
        hw.add_control_port(ctrl_port)

    # Paper (4.3)-(4.5): DG_k heat scalar, RT_k heat flux,
    # CG_(k+1) interface traces. Our existing wave grad-grad PFEM
    # requires a continuous velocity p: CG_(k+1) rather than the
    # article's div-div DG_k. Wave stress q is now H(div)-conforming
    # RT_k, giving a genuine normal trace.
    # k=0 is the order used in the reference numerical experiments.
    k = 0
    FEMs = [
        S.FEM("T", k, "DG"),
        S.FEM("Heat flux", k, "RT"),
        S.FEM("Interface Heat", k + 1, "CG"),
        S.FEM("p", k + 1, "CG"),
        S.FEM("q", k, "RT"),
        S.FEM("Interface Wave", k + 1, "CG"),
        S.FEM("Boundary", k + 1, "CG"),
    ]
    
    # Add them to the dphs
    for FEM in FEMs:
        hw.add_FEM(FEM)
        # Earlier GetFEM versions could silently fall back to CG
        # on an unavailable intrinsic RT element. Never claim an
        # RT-based model if that fallback occurred.
        if FEM.get_name() in ("Heat flux", "q") and FEM.get_type() != "RT":
            raise RuntimeError(
                "RT finite elements are required for the heat-wave example"
            )

    # Define the pHs via `Brick` == non-zero block matrices == variational terms
    # Since we use co-energy formulation, constitutive relations are already taken into
    # account in the mass matrices M_q and M_p
    bricks = [
        # === Heat: div-div
        S.Brick("M_T", "T*Test_T", [heat_region], dt=True, position="flow"),
        S.Brick("M_Q", "e_Q.Test_e_Q", [heat_region], position="flow"),
        S.Brick("M_Y_T", "Y_T*Test_Y_T", [10], position="flow"),
        
        S.Brick("D_T", "-Div(e_Q)*Test_T", [heat_region], position="effort"),
        S.Brick("D_T^T", "T*Div(Test_e_Q)", [heat_region], position="effort"),
        # With e_Q=-grad(T) and the outward HEAT normal, T=U_T
        # enters with a minus in the mixed boundary weak form,
        # for either assignment of heat_region.
        S.Brick(
            "B_T", "-U_T*Test_e_Q.Normal",
            [10], position="effort",
        ),
        S.Brick("B_T^T", "e_Q.Normal*Test_Y_T", [10], position="effort"),

        # === Wave: grad-grad
        S.Brick("M_p", "p*Test_p", [wave_region], dt=True, position="flow"),
        S.Brick("M_q", "q.Test_q", [wave_region], dt=True, position="flow"),
        S.Brick("M_Y_w", "Y_w*Test_Y_w", [10], position="flow"),
        
        S.Brick("D_w", "-q.Grad(Test_p)", [wave_region], position="effort"),
        S.Brick("-D_w^T", "Grad(p).Test_q", [wave_region], position="effort"),
        S.Brick("B_w", "U_w*Test_p", [10], position="effort"),
        S.Brick("B_w^T", "p*Test_Y_w", [10], position="effort"),
    ]
    # Outer homogeneous Dirichlet, consistently with (2.1):
    # - heat exterior: prescribe temperature, observe normal flux;
    # - wave exterior: prescribe velocity, normal traction is the
    #   Lagrange multiplier, following examples/wave.py.
    if not outer_is_wave:
        bricks.append(S.Brick(
            "M_Y_bnd", "Y_bnd*Test_Y_bnd", [20], position="flow"
        ))
        bricks.append(S.Brick(
            "B_bnd", "-U_bnd*Test_e_Q.Normal", [20], position="effort"
        ))
        bricks.append(S.Brick(
            "B_bnd^T", "e_Q.Normal*Test_Y_bnd", [20], position="effort"
        ))
    else:
        bricks.append(S.Brick(
            "M_Y_bnd", "U_bnd*Test_Y_bnd", [20], position="flow"
        ))
        bricks.append(S.Brick(
            "B_bnd", "Y_bnd*Test_p", [20], position="effort"
        ))
        bricks.append(S.Brick(
            "B_bnd^T", "-p*Test_Y_bnd", [20], position="effort"
        ))
    for brick in bricks:
        hw.add_brick(brick)
    
    # Set the controls
    # Physical continuity (2.2): T=p, J_Q.n_heat=q.n_wave, with
    # each Normal outward from its OWN subdomain. The two interface
    # controls impose U_T=Y_w and U_w=Y_T in both region assignments.
    # For the mixed/grad-grad formulation above, the signed interface
    # power is -U_T*Y_T + U_w*Y_w = 0. No disk/annulus sign factor.
    hw.set_control("Interface Heat", "Y_w")
    hw.set_control("Interface Wave", "Y_T")
    # === Dirichlet boundary condition
    hw.set_control("Boundary", "0.")

    # All initial fields use the actual FE L2 projection. The original
    # Gaussian is the physical-reference case. The optional gallery
    # excitation is a boundary-localized angular wave packet, NOT a
    # claimed whispering-gallery eigenfunction. Both preserve the same
    # coupled PDE and the exactly homogeneous boundary condition.
    if initial_wave_profile not in ("gaussian", "gallery_packet"):
        raise ValueError("initial_wave_profile must be gaussian or gallery_packet")
    if initial_wave_profile == "gallery_packet":
        if wave_region != 2:
            raise ValueError("gallery_packet requires wave_region=2 (annulus)")
        # r in (0.6, 1.0), theta around the annulus. The continuous
        # packet vanishes at r=0.6 and r=1, but unconstrained L2
        # projection may reintroduce nonzero essential boundary DOFs.
        # The same qualified discrete clamp below still applies.
        radius = "np.sqrt(x*x+y*y)"
        p_initial = (
            f"50.*({radius}-0.6)*(1.-{radius})*"
            f"np.exp(-120.*({radius}-0.86)**2)*"
            "np.cos(14.*np.arctan2(y,x))"
        )
        T_initial = "0."
    else:
        T_initial = "5.*np.exp(-25*((x-0.6)*(x-0.6)+y*y))"
        p_initial = T_initial
    hw.set_initial_value("T", T_initial)
    hw.set_initial_value("p", p_initial)
    hw.set_initial_value("q", "[0.,0.]")

    # For a wave annulus, the Gaussian projected in CG1 has nonzero
    # exterior-boundary coefficients, while homogeneous Dirichlet
    # requires p=0 there. The strict algebraic DAE initializer must
    # REJECT that initial field, not silently alter its dynamic state.
    #
    # Both remedies were previously native-qualified on this exact
    # Concentric Heat-Wave system: (1) clamp only boundary CG1 DOFs,
    # or (2) mass-orthogonal constrained L2 correction. The default
    # is the minimal transparent clamp; both keep T and q unchanged.
    # When heat is outside, there is no exterior wave Dirichlet trace.
    if initial_wave_boundary not in ("clamp", "constrained_l2", "unconstrained"):
        raise ValueError("initial_wave_boundary must be clamp, constrained_l2 "
                         "or unconstrained")
    hw.heat_wave_initialization_report = {
        "method": "not_required_heat_outer" if not outer_is_wave
                  else "unconstrained_historical_diagnostic",
        "boundary_compatibility_enforced": False,
    }
    hw.heat_wave_initialization_profile = initial_wave_profile
    if outer_is_wave and initial_wave_boundary != "unconstrained":
        original_T = np.asarray(hw.gf_model.variable("T")).reshape(-1).copy()
        original_q = np.asarray(hw.gf_model.variable("q")).reshape(-1).copy()
        report = enforce_zero_dirichlet_initial_values(
            hw, "p", volume_region=wave_region, boundary_region=20,
            method=initial_wave_boundary, factor_solver="mumps",
        )
        if (not np.array_equal(
                np.asarray(hw.gf_model.variable("T")).reshape(-1), original_T
            ) or not np.array_equal(
                np.asarray(hw.gf_model.variable("q")).reshape(-1), original_q
            )):
            raise RuntimeError(
                "Wave boundary compatibility changed the heat/stress initial data"
            )
        report["boundary_compatibility_enforced"] = True
        hw.heat_wave_initialization_report = report

    ## Solve in time
    # CN is the second-order time integrator used for the paper's
    # energy-decay results; fixed accepted steps make the model's
    # dissipation work independently checkable step by step.
    # No explicit solver => SCRIMP's qualified MUMPS default.
    # Other PETSc LU backends remain an opt-in user choice.
    scheme = {}
    if factor_solver is not None:
        scheme["pc_factor_mat_solver_type"] = factor_solver
    hw.set_time_scheme(
        # Use the production constraint-consistent initializer for the
        # index-2 coupled pH-DAE; do not hide an inconsistent Gaussian.
        init_step_ts_type="algebraic",
        ts_type="cn",
        ts_adapt_type="none",
        t_f=float(final_time),
        dt=float(time_step),
        dt_save=float(save_step),
                       ksp_type="preonly",
                       pc_type="lu",
                       **scheme,
                       )

    # Solve
    hw.solve()

    ## Post-processing
    ## Set Hamiltonian's name
    hw.hamiltonian.set_name("Hamiltonian")
    # Define each Hamiltonian Term
    terms = [
        S.Term("Lyapunov heat", "0.5*T*T", [heat_region]),
        S.Term("Kinetic energy", "0.5*p*p", [wave_region]),
        S.Term("Potential energy", "0.5*q.q", [wave_region]),
    ]
    # Add them to the Hamiltonian
    for term in terms:
        hw.hamiltonian.add_term(term)

    # Plot the Hamiltonian and save the output
    hw.plot_Hamiltonian(save_figure=True, filename="Hamiltonian_Heat"+str(heat_region)+"_Wave"+str(wave_region)+"_2D.png")
    
    # A true log-log diagnostic must omit t=0, rather than plotting
    # zero on a logarithmic axis. Keep this separate from the energy
    # identity plot used for the physical qualification.
    import matplotlib.pyplot as plt
    t = np.asarray(hw.solution["t"], dtype=float)
    Hamiltonian = np.asarray(hw.get_Hamiltonian(), dtype=float)
    fig, ax = plt.subplots(figsize=(8, 5))
    positive = (t > 0.) & (Hamiltonian > 0.)
    ax.loglog(t[positive], Hamiltonian[positive])
    ax.grid(axis="both")
    ax.set_xlabel("time t")
    ax.set_ylabel("Hamiltonian")
    ax.set_title("Evolution of the Hamiltonian (log-log)")
    # The example itself saves this long-time diagnostic, with t=0
    # excluded because log(0) is undefined. Never fit an exponent
    # unless a separate physical convergence/decay study validates it.
    from pathlib import Path
    png = Path.cwd()/"outputs"/"png"
    png.mkdir(parents=True,exist_ok=True)
    fig.tight_layout()
    fig.savefig(png/("Hamiltonian_Heat"+str(heat_region)+
                     "_Wave"+str(wave_region)+"_2D_loglog.png"),dpi=170)
    fig.savefig(png/("Hamiltonian_Heat"+str(heat_region)+
                     "_Wave"+str(wave_region)+"_2D_loglog.svg"))
    plt.show()
    plt.close(fig)
    
    hw.export_to_pv("p")
    hw.export_to_pv("T")
    if long_showcase:
        if wave_region != 2:
            raise ValueError("Long annular showcase requires wave_region=2")
        # Original genuine RT stress field, also browseable in ParaView.
        hw.export_to_pv("q")
        hw.heat_wave_showcase_report = render_heat_wave_long_showcase(
            t, Hamiltonian,
            profile=initial_wave_profile, wave_region=wave_region,
        )

    return hw

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Canonical native Heat–Wave 2D pH simulation, "
                    "with optional long-time annular wave visualizations"
    )
    # Historical two-integer invocation stays supported.
    parser.add_argument("heat_positional",nargs="?",type=int)
    parser.add_argument("wave_positional",nargs="?",type=int)
    parser.add_argument("--heat-region",type=int,default=1)
    parser.add_argument("--wave-region",type=int,default=2)
    parser.add_argument("--final-time",type=float,default=2.0)
    parser.add_argument("--time-step",type=float,default=0.001)
    parser.add_argument("--save-step",type=float,default=0.05)
    parser.add_argument("--factor-solver",choices=("mumps","superlu_dist"),
                        default=None)
    parser.add_argument("--initial-wave-boundary",
                        choices=("clamp","constrained_l2","unconstrained"),
                        default="clamp")
    parser.add_argument("--initial-wave-profile",
                        choices=("gaussian","gallery_packet"),
                        default="gaussian")
    parser.add_argument("--long-showcase",action="store_true")
    args=parser.parse_args()
    if (args.heat_positional is None)!=(args.wave_positional is None):
        parser.error("Both historical heat/wave positional regions are required")
    heat_region = (args.heat_positional if args.heat_positional is not None
                   else args.heat_region)
    wave_region = (args.wave_positional if args.wave_positional is not None
                   else args.wave_region)
    if (not np.isfinite(args.final_time) or args.final_time<=0
            or not np.isfinite(args.time_step) or args.time_step<=0
            or not np.isfinite(args.save_step) or args.save_step<=0):
        parser.error("Times and steps must be positive finite numbers")
    hw = heat_wave_eq(
        heat_region=heat_region, wave_region=wave_region,
        final_time=args.final_time, time_step=args.time_step,
        save_step=args.save_step, factor_solver=args.factor_solver,
        initial_wave_boundary=args.initial_wave_boundary,
        initial_wave_profile=args.initial_wave_profile,
        long_showcase=args.long_showcase,
    )
    
