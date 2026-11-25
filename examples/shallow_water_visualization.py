# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Real GetFEM/VTK scientific pictures for canonical inviscid Saint-Venant.

Imported and RUN ONLY when examples/shallow_water.py requests --showcase.
All spatial values are read directly from the h (CG2) and p (DG1)
native GetFEM VTK exports, with no FE indexing assumptions, no synthetic
water heights and no interpolation between times. Velocity is p / rho,
NOT raw p. A cellwise triangular average serves only as a plotting
convention; original VTK arrays are retained for ParaView/verification.

Native GetFEM + BDF2 physical acceptance and visualization acceptance
are separate: no artificial exact Hamiltonian conservation claim.
"""
import base64
import io
import json
from pathlib import Path
import struct
from xml.etree import ElementTree as ET

import numpy as np


_TYPE = {
    "Float32": "<f4", "Float64": "<f8",
    "Int32": "<i4", "Int64": "<i8",
    "UInt8": "u1", "UInt32": "<u4", "UInt64": "<u8",
}


def vtk_array(node):
    """Parse genuine GetFEM ASCII / inline-binary VTK numeric DataArray."""
    if node is None:
        raise ValueError("Missing original GetFEM VTU DataArray")
    dtype_name = node.get("type", "")
    if dtype_name not in _TYPE:
        raise ValueError("Unexpected native VTK array type: "+dtype_name)
    dtype = np.dtype(_TYPE[dtype_name])
    fmt = node.get("format", "ascii")
    if fmt == "ascii":
        data = np.fromstring(node.text or "", dtype=dtype, sep=" ")
    elif fmt == "binary":
        payload = base64.b64decode(
            "".join((node.text or "").split()), validate=True
        )
        if len(payload) < 4:
            raise ValueError("Truncated native GetFEM VTU data")
        nbytes = struct.unpack("<I", payload[:4])[0]
        if nbytes != len(payload)-4 or nbytes % dtype.itemsize:
            raise ValueError("Invalid binary VTK byte count")
        data = np.frombuffer(payload, dtype=dtype, offset=4).copy()
    else:
        raise ValueError("Unsupported GetFEM VTK encoding: "+fmt)
    if not np.all(np.isfinite(data)):
        raise ValueError("Nonfinite actual GetFEM field/mesh value")
    return data


def native_triangular_field(filename, field_name):
    """Return original GetFEM triangle points, cells, h/p field values."""
    filename = Path(filename)
    piece = ET.parse(filename).getroot().find(".//Piece")
    if piece is None:
        raise ValueError("Missing native GetFEM VTU Piece")
    npoints = int(piece.attrib.get("NumberOfPoints", "0"))
    ncells = int(piece.attrib.get("NumberOfCells", "0"))
    if npoints<3 or ncells<1:
        raise ValueError("Empty original GetFEM triangle mesh")
    points = vtk_array(piece.find("./Points/DataArray"))
    if points.shape != (3*npoints,):
        raise ValueError("GetFEM VTK coordinates not 3D embedded")
    points = points.reshape(npoints,3)
    if np.max(abs(points[:,2]))>1e-7:
        raise ValueError("Expected planar original shallow-water rectangle")

    con=vtk_array(piece.find("./Cells/DataArray[@Name='connectivity']"))
    off=vtk_array(piece.find("./Cells/DataArray[@Name='offsets']"))
    typ=vtk_array(piece.find("./Cells/DataArray[@Name='types']"))
    if (con.shape!=(3*ncells,) or off.shape!=(ncells,)
            or typ.shape!=(ncells,)
            or not np.array_equal(off,3*np.arange(1,ncells+1))
            or not np.all(typ==5)):
        raise ValueError("Expected actual linear VTK triangles from GetFEM")
    triangles=np.asarray(con.reshape(ncells,3),dtype=np.int64)
    if np.min(triangles)<0 or np.max(triangles)>=npoints:
        raise ValueError("Invalid original GetFEM VTK connectivity")

    node=piece.find(f"./PointData/DataArray[@Name='{field_name}']")
    components = int(node.get("NumberOfComponents","1")) if node is not None else 0
    if field_name=="h" and components!=1:
        raise ValueError("Missing scalar free surface h in GetFEM VTU")
    if field_name=="p" and components not in (2,3):
        raise ValueError("Missing vector original momentum p in GetFEM VTU")
    data=vtk_array(node)
    if field_name=="h":
        if data.shape!=(npoints,):
            raise ValueError("Invalid native h CG2 field size")
    else:
        if data.size!=npoints*components:
            raise ValueError("Invalid native DG1 momentum vector size")
        data=data.reshape(npoints,components)[:,:2]
    return points[:,:2].astype(float),triangles,data.astype(float)


def true_pvd_files(root, name, true_times):
    """Demand exactly one genuine native VTK snapshot for each saved FE t."""
    pvd = root/"pv"/name/(name+".pvd")
    if not pvd.is_file():
        raise FileNotFoundError("Missing original GetFEM PVD: "+str(pvd))
    sets=ET.parse(pvd).getroot().findall(".//DataSet")
    if len(sets)!=len(true_times):
        raise ValueError("GetFEM PVD and true FE solution lengths disagree")
    result=[]
    for index,node in enumerate(sets):
        timestamp=float(node.attrib["timestep"])
        if abs(timestamp-true_times[index])>1e-8:
            raise ValueError("GetFEM PVD time differs from original FE time")
        relative=Path(node.attrib["file"])
        if (relative.name!=str(relative)
                or relative.suffix!=".vtu"):
            raise ValueError("Unsafe / non-original GetFEM VTK filename")
        filename=pvd.parent/relative
        if not filename.is_file() or filename.stat().st_size<100:
            raise ValueError("Original GetFEM VTK field missing: "+str(filename))
        result.append(filename)
    return result


def render_native_shallow_water(system, folder, *, rho=1000.,
                               max_gif_frames=36):
    """Render true 2-D h and p/rho from exported VTK, plus physical budgets.

    Call on the actual canonical DPHS instance AFTER successful
    nonlinear_affine/BDF2 simulation, H registration, and native
    export_to_pv('h')/export_to_pv('p'). Only MPI1 is supported for
    this visual protocol; MPI2 physical parity is a separate audit.
    """
    from matplotlib import pyplot as plt
    from matplotlib import cm, colors
    from matplotlib.collections import PolyCollection
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from PIL import Image

    folder=Path(folder).resolve()
    folder.mkdir(parents=True,exist_ok=True)
    figures=folder/"figures"
    figures.mkdir(exist_ok=True)
    t=np.asarray(system.solution["t"],dtype=float).reshape(-1)
    H=np.asarray(system.get_Hamiltonian(),dtype=float).reshape(-1)
    mass=np.asarray(system.get_quantity("h",region=1),dtype=float).reshape(-1)
    kinetic=np.asarray(
        system.get_quantity("0.5*h*p.p/rho",region=1),dtype=float
    ).reshape(-1)
    potential=np.asarray(
        system.get_quantity("0.5*rho*g*h*h",region=1),dtype=float
    ).reshape(-1)
    contributions=[]
    for i in range(4):
        contributions.append(
            np.asarray(system.get_quantity(
                f"U_{i}*Y_{i}",region=10+i
            ),dtype=float).reshape(-1)
        )
    P=np.sum(contributions,axis=0)
    if (len(t)<5 or H.shape!=t.shape or
            any(arr.shape!=t.shape for arr in (mass,kinetic,potential,P))
            or not np.all(np.isfinite(t)) or not np.all(np.diff(t)>0)
            or not all(np.all(np.isfinite(a)) for a in
                       (H,mass,kinetic,potential,P))
            or np.any(H<=0) or np.any(mass<=0)
            or np.max(abs(H-kinetic-potential))>1e-6):
        raise ValueError("Invalid authentic GetFEM BDF2 physical time traces")
    boundary_work=np.zeros_like(t)
    boundary_work[1:]=np.cumsum(
        0.5*(P[1:]+P[:-1])*np.diff(t)
    )
    energy_defect=H-H[0]-boundary_work
    h_files=true_pvd_files(folder,"h",t)
    p_files=true_pvd_files(folder,"p",t)

    # Inspect each actual snapshot ONCE, so fixed ranges represent the
    # global native FE extrema, not arbitrary per-frame normalization.
    selected=np.unique(np.linspace(
        0,len(t)-1,min(len(t),max_gif_frames),dtype=int
    ))
    h_frames=[]
    p_frames=[]
    h_geometry=None
    p_geometry=None
    hmin,hmax=float("inf"),-float("inf")
    speedmax=0.
    for idx in selected:
        hx,ht,h=native_triangular_field(h_files[int(idx)],"h")
        px,pt,p=native_triangular_field(p_files[int(idx)],"p")
        if h_geometry is None:
            h_geometry=(hx,ht)
            p_geometry=(px,pt)
        elif (not np.array_equal(h_geometry[0],hx)
              or not np.array_equal(h_geometry[1],ht)
              or not np.array_equal(p_geometry[0],px)
              or not np.array_equal(p_geometry[1],pt)):
            raise ValueError("Original GetFEM VTK geometry changed across time")
        hmin=min(hmin,float(np.min(h)))
        hmax=max(hmax,float(np.max(h)))
        if np.min(h)<=0:
            raise ValueError("Nonphysical nonpositive native water height")
        speedmax=max(speedmax,float(np.max(np.linalg.norm(p,axis=1))/rho))
        h_frames.append(h)
        p_frames.append(p)
    if hmax<=hmin or speedmax<=0:
        raise ValueError("Native dam-break wave has no nontrivial h or velocity")

    # Original GetFEM triangular cells directly (CG2 h and DG1 p
    # use independent VTK point sets). Cell averages are ONLY
    # per-triangle graphics, not a substitute FE mass or energy.
    hxy,htri=h_geometry
    pxy,ptri=p_geometry
    h_cells=hxy[htri]
    p_centres=pxy[ptri].mean(axis=1)
    p_step=max(1,len(ptri)//160)
    hnorm=colors.Normalize(vmin=hmin,vmax=hmax)
    vnorm=colors.Normalize(vmin=0.,vmax=speedmax)
    h_cmap=cm.get_cmap("viridis")
    v_cmap=cm.get_cmap("plasma")

    def draw_frame(k):
        height=h_frames[k]
        momentum=p_frames[k]
        height_cells=height[htri].mean(axis=1)
        vel=momentum[ptri].mean(axis=1)/rho
        speeds=np.linalg.norm(vel,axis=1)
        fig=plt.figure(figsize=(14.2,4.9),constrained_layout=True)
        grid=fig.add_gridspec(1,3,width_ratios=(1.1,1.5,1.2))
        ax_surface=fig.add_subplot(grid[0,0],projection="3d")
        ax_height=fig.add_subplot(grid[0,1])
        ax_velocity=fig.add_subplot(grid[0,2])
        # Original VTK triangle vertices at their actual nodal h.
        # 3D geometry is a *visual elevation* of the 2-D FE surface,
        # NOT a separate 3-D shallow-water simulation.
        xyz=np.stack((
            h_cells[:,:,0],h_cells[:,:,1],height[htri]
        ),axis=2)
        surface=Poly3DCollection(
            xyz,linewidth=0.,edgecolor="none"
        )
        surface.set_facecolor(h_cmap(hnorm(height_cells)))
        ax_surface.add_collection3d(surface)
        ax_surface.set_xlim(0,2)
        ax_surface.set_ylim(0,.5)
        ax_surface.set_zlim(max(0,hmin-.15),hmax+.15)
        ax_surface.set_box_aspect((2.,.5,.65),zoom=.9)
        ax_surface.view_init(elev=26,azim=-68)
        ax_surface.set(xlabel="x",ylabel="y",zlabel="h")
        ax_surface.set_title("Surface libre 3D")

        collection=PolyCollection(
            h_cells, array=height_cells, cmap=h_cmap, norm=hnorm,
            edgecolors="none",linewidths=0,
        )
        ax_height.add_collection(collection)
        ax_height.set(aspect="equal",xlim=(0.,2.),ylim=(0.,.5),
                      xlabel="x",ylabel="y",
                      title="Hauteur d'eau h (champ GetFEM)")
        fig.colorbar(cm.ScalarMappable(norm=hnorm,cmap=h_cmap),
                     ax=ax_height,orientation="horizontal",
                     shrink=.76,label="h")

        sample=np.arange(0,len(ptri),p_step)
        quiver=ax_velocity.quiver(
            p_centres[sample,0],p_centres[sample,1],
            vel[sample,0],vel[sample,1],
            speeds[sample],cmap=v_cmap,norm=vnorm,
            angles="xy",scale_units="xy",
            scale=max(1e-12,3.*speedmax),
            pivot="middle",width=.0035,
        )
        ax_velocity.set(
            aspect="equal",xlim=(0.,2.),ylim=(0.,.5),
            xlabel="x",ylabel="y",
            title="Vitesse v = p/ρ (DG1 GetFEM)"
        )
        fig.colorbar(cm.ScalarMappable(norm=vnorm,cmap=v_cmap),
                     ax=ax_velocity,orientation="horizontal",
                     shrink=.76,label="||v||")
        fig.suptitle(
            "SCRIMP v1.3.0 · Saint–Venant 2D non visqueux · "
            f"t={t[selected[k]]:.4f} · vrais champs FE",
            fontsize=12,
        )
        return fig

    # Scientific plots from the exact FE integral histories, NOT from
    # the plotting-only cell averages or original VTK Float32 values.
    fig,axes=plt.subplots(2,1,figsize=(9,7),sharex=True)
    axes[0].plot(t,H,label="H mécanique GetFEM")
    axes[0].plot(t,kinetic,label="Énergie cinétique")
    axes[0].plot(t,potential,label="Énergie potentielle")
    axes[0].set_ylabel("Énergies");axes[0].legend()
    axes[0].grid(alpha=.3)
    axes[1].plot(t,energy_defect,label="H(t)−H(0)−W_bord(t)")
    axes[1].plot(t,(mass-mass[0])/mass[0],
                 label="Dérive relative de masse")
    axes[1].set(xlabel="Temps",ylabel="Défauts",title=
                "BDF2 : bilan mesuré, sans conservation exacte postulée")
    axes[1].legend();axes[1].grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(figures/"saint_venant_H_masse_bilan.png",dpi=170)
    fig.savefig(figures/"saint_venant_H_masse_bilan.svg")
    plt.close(fig)

    snapshots=sorted(set((0,len(selected)//3,
                          2*len(selected)//3,len(selected)-1)))
    for k in snapshots:
        fig=draw_frame(k)
        fig.savefig(
            figures/f"saint_venant_champs_t{t[selected[k]]:.4f}.png",
            dpi=160,
        )
        plt.close(fig)
    frames=[]
    for k in range(len(selected)):
        fig=draw_frame(k)
        buffer=io.BytesIO()
        fig.savefig(buffer,format="png",dpi=90)
        plt.close(fig)
        buffer.seek(0)
        with Image.open(buffer) as image:
            frames.append(image.convert("RGB").copy())
    frames[0].save(
        figures/"saint_venant_surface_vitesse.gif",
        save_all=True,append_images=frames[1:],
        duration=160,loop=0,optimize=False,
    )

    np.savez_compressed(
        folder/"physical_FE_traces.npz",
        t=t,H=H,mass=mass,kinetic=kinetic,potential=potential,
        boundary_power=P, boundary_work=boundary_work,
        numerical_BDF2_energy_defect=energy_defect,
    )
    provenance={
        "schema":"v130-original-inviscid-shallow-water-native-visuals-v1",
        "status":"NATIVE_GETFEM_PVD_VISUALS_NOT_SPATIAL_CONVERGENCE",
        "model":"examples/shallow_water.py",
        "velocity_definition":"v=p/rho",
        "rho":float(rho),
        "FEM_h":"CG2","FEM_p":"DG1",
        "PDE_modified":False,
        "initialization":"nonlinear_affine (explicit showcase opt-in)",
        "integrator":"BDF2, explicit lagged gyroscopic G",
        "H_exact_conservation_claimed":False,
        "MPI2_visuals_qualified":False,
        "viscous_rotating_research_included":False,
        "source_of_spatial_values":"original GetFEM h/p native PVD and VTU",
        "plotting_rule":"mean of each actual VTK triangle's three native values",
        "no_time_interpolation":True,
        "fixed_height_color_range":[float(hmin),float(hmax)],
        "fixed_velocity_color_range":[0.,float(speedmax)],
        "true_FE_states":int(len(t)),"gif_original_FE_frames":int(len(selected)),
        "original_selected_times":t[selected].tolist(),
        "mesh_h_native_triangles":int(len(htri)),
        "mesh_p_native_triangles":int(len(ptri)),
        "Hamiltonian_initial":float(H[0]),
        "Hamiltonian_final":float(H[-1]),
        "max_relative_mass_drift":float(
            np.max(abs(mass-mass[0]))/mass[0]
        ),
        "final_BDF2_energy_defect":float(energy_defect[-1]),
        "max_original_boundary_power":float(np.max(abs(P))),
    }
    (folder/"provenance.json").write_text(
        json.dumps(provenance,indent=2)+"\n",encoding="utf-8",
    )
    return provenance
