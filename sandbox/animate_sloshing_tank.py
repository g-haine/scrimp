# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Render a moving-tank and free-surface MP4 from SCRIMP solution snapshots.

Usage, from the repository root::

    python sandbox/animate_sloshing_tank.py --output outputs/sloshing_tank.mp4

The numerical model remains the experimental v1.2.0 sandbox model. All
fields are obtained from GetFEM/SCRIMP; only intermediate movie frames
are linearly interpolated in time. Display exaggerations are labeled.
"""

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class MovieData:
    """Sampled SCRIMP solution in physical units before display scaling."""

    times: np.ndarray               # (nt,)
    x: np.ndarray                   # (nx,), tank-relative coordinate
    depths: np.ndarray              # (nt, nx), physical water depth
    positions: np.ndarray           # (nt,), physical tank displacement
    velocities: np.ndarray          # (nt,), physical tank velocity
    energies: np.ndarray            # (nt,), total Hamiltonian
    input_work: np.ndarray          # (nt,), integrated applied power
    viscous_loss: np.ndarray        # (nt,), integrated viscous dissipation
    balance_defect: np.ndarray      # (nt,), dH - (work - dissipation)
    input_power: np.ndarray         # (nt,), applied force * tank velocity
    viscous_power: np.ndarray       # (nt,), nonnegative dissipation rate
    time_scheme: str = "cn"


def _real_array(data, name):
    """Read an effectively real quantity without silently dropping imaginary data."""
    values = np.asarray(data)
    if np.iscomplexobj(values):
        scale = max(1.0, float(np.max(np.abs(values.real))))
        if np.max(np.abs(values.imag)) > 1e-9 * scale:
            raise ValueError(f"{name} contains significant imaginary values")
        values = values.real
    values = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} contains non-finite values")
    return values


def validate_data(data):
    """Fail early on malformed or unphysical data instead of drawing it."""
    t = np.asarray(data.times)
    x = np.asarray(data.x)
    if t.ndim != 1 or len(t) < 2 or not np.all(np.diff(t) > 0):
        raise ValueError("At least two strictly increasing saved times are needed")
    if x.ndim != 1 or len(x) < 2 or not np.all(np.diff(x) > 0):
        raise ValueError("x must contain increasing tank-relative coordinates")
    if data.depths.shape != (len(t), len(x)):
        raise ValueError("depths must have shape (len(times), len(x))")
    for name in ("positions", "velocities", "energies",
                 "input_work", "viscous_loss", "balance_defect",
                 "input_power", "viscous_power"):
        if np.asarray(getattr(data, name)).shape != (len(t),):
            raise ValueError(f"{name} must have shape (len(times),)")
    for name in ("times", "x", "depths", "positions", "velocities",
                 "energies", "input_work", "viscous_loss", "balance_defect",
                 "input_power", "viscous_power"):
        if not np.all(np.isfinite(getattr(data, name))):
            raise ValueError(f"{name} contains non-finite values")
    if np.min(data.depths) <= 0:
        raise ValueError("The computed water depth is non-positive")
    residual = (data.energies - data.energies[0]) - (data.input_work - data.viscous_loss)
    if not np.allclose(residual, data.balance_defect, atol=1e-9, rtol=1e-6):
        raise ValueError("Energy/power-balance time series are inconsistent")


def extract_movie_data(tank, samples=250):
    """Interpolate actual saved finite-element area coefficients on a 1D grid.

    Requires a completed sandbox.sloshing_tank.sloshing_tank run.
    GetFEM's documented compute_interpolate_on provides values at
    arbitrary points, including for the current quadratic CG area.
    """
    import getfem as gf

    if samples < 3:
        raise ValueError("samples must be >= 3")
    geometry = getattr(tank, "sloshing_visualization", None)
    if geometry is None:
        raise ValueError("Missing geometry metadata: use the updated sandbox/sloshing_tank.py")
    length = float(geometry["length"])
    width = float(geometry["width"])
    if length <= 0 or width <= 0:
        raise ValueError("Tank length and width must be positive")

    times = _real_array(tank.solution["t"], "times").ravel()
    x = np.linspace(0.0, length, int(samples))
    points = x.reshape(1, -1)
    mesh_fem = tank.ports["a_h"].get_fem()
    area_snapshots = tank.get_solution("a_h")
    if len(area_snapshots) != len(times):
        raise ValueError("Inconsistent area snapshots and recorded times")

    depths = []
    for n, area_dofs in enumerate(area_snapshots):
        # Finite-element area a_h = width * h, NOT a prescribed sine wave.
        coeffs = _real_array(area_dofs, f"a_h coefficients at step {n}").ravel()
        interpolated = gf.compute_interpolate_on(mesh_fem, coeffs, points)
        area = _real_array(interpolated, f"a_h at step {n}").ravel()
        if len(area) != len(x):
            raise ValueError(f"GetFEM returned {len(area)} points; expected {len(x)}")
        depths.append(area / width)

    diagnostic = getattr(tank, "sloshing_power_balance", None)
    if diagnostic is None:
        raise ValueError("Missing power balance: use the updated sandbox/sloshing_tank.py")
    if not np.allclose(_real_array(diagnostic["times"], "balance times"), times):
        raise ValueError("Power balance does not correspond to saved times")

    data = MovieData(
        times=times,
        x=x,
        depths=np.asarray(depths),
        positions=_real_array(tank.get_solution("pos"), "pos").ravel(),
        velocities=_real_array(tank.get_solution("v"), "v").ravel(),
        energies=_real_array(diagnostic["hamiltonian"], "Hamiltonian").ravel(),
        input_work=_real_array(diagnostic["input_work"], "input work").ravel(),
        viscous_loss=_real_array(diagnostic["viscous_loss"], "viscous loss").ravel(),
        balance_defect=_real_array(diagnostic["defect"], "balance defect").ravel(),
        input_power=_real_array(diagnostic["input_power"], "input power").ravel(),
        viscous_power=_real_array(diagnostic["viscous_dissipation_power"], "viscous power").ravel(),
        time_scheme=str(diagnostic["time_scheme"]),
    )
    validate_data(data)
    return data



def save_balance_csv(data, output):
    """Export computed time histories (physical units, no movie exaggeration)."""
    validate_data(data)
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("time_s", "hamiltonian", "delta_h", "input_power",
                         "viscous_power", "input_work", "viscous_loss",
                         "expected_delta_h", "balance_defect", "time_scheme"))
        for i, t in enumerate(data.times):
            writer.writerow((t, data.energies[i],
                             data.energies[i]-data.energies[0],
                             data.input_power[i], data.viscous_power[i],
                             data.input_work[i], data.viscous_loss[i],
                             data.input_work[i]-data.viscous_loss[i],
                             data.balance_defect[i], data.time_scheme))
    return out


def interpolate_frame(data, t):
    """Linear interpolation between two *saved* numerical states for display."""
    times = data.times
    if t <= times[0]:
        return data.depths[0], float(data.positions[0]), float(data.velocities[0]), float(data.energies[0])
    if t >= times[-1]:
        return data.depths[-1], float(data.positions[-1]), float(data.velocities[-1]), float(data.energies[-1])
    hi = int(np.searchsorted(times, t))
    lo = hi - 1
    alpha = (t - times[lo]) / (times[hi] - times[lo])
    return (
        data.depths[lo] + alpha * (data.depths[hi] - data.depths[lo]),
        float(data.positions[lo] + alpha * (data.positions[hi] - data.positions[lo])),
        float(data.velocities[lo] + alpha * (data.velocities[hi] - data.velocities[lo])),
        float(data.energies[lo] + alpha * (data.energies[hi] - data.energies[lo])),
    )


def render_movie(data, output, fps=25, seconds=10.0, surface_scale=8.0,
                 motion_scale="auto", dpi=120):
    """Write a headless MP4 with water, moving tank and diagnostic plots.

    The surface perturbation (h - mean(h)) is display-scaled, not h itself.
    Tank motion may also be display-scaled. Numeric readouts and the two
    time-series plots always show *unscaled physical* values.

    The power balance is a *measured numerical defect*, not a
    guarantee of exact energy preservation for nonlinear CN integration.
    """
    import matplotlib
    matplotlib.use("Agg")  # Save directly to MP4, without Qt/Wayland.
    from matplotlib import pyplot as plt, animation
    from matplotlib.patches import Polygon

    validate_data(data)
    if fps < 1 or seconds <= 0 or dpi < 40 or surface_scale <= 0:
        raise ValueError("fps, seconds, dpi and surface_scale must be positive")
    if not animation.writers.is_available("ffmpeg"):
        raise RuntimeError("FFmpeg not found: conda install -c conda-forge ffmpeg")

    output = Path(output)
    if output.suffix.lower() != ".mp4":
        raise ValueError("The output filename must end in .mp4")
    output.parent.mkdir(parents=True, exist_ok=True)

    length = float(data.x[-1] - data.x[0])
    excursion = np.max(np.abs(data.positions - data.positions[0]))
    if motion_scale == "auto":
        move_factor = max(1.0, 0.16 * length / excursion) if excursion > 1e-15 else 1.0
    else:
        move_factor = float(motion_scale)
        if not np.isfinite(move_factor) or move_factor <= 0:
            raise ValueError("motion_scale must be 'auto' or a positive number")

    # Fixed physical reference level: scale surface deviations only.
    reference_height = float(np.sum(
        0.5 * (data.depths[0][1:] + data.depths[0][:-1]) * np.diff(data.x)
    ) / length)
    displayed_depths = reference_height + surface_scale * (data.depths - reference_height)
    if np.min(displayed_depths) <= 0:
        raise ValueError("surface_scale is too high: displayed water would cross the tank floor")

    offset = (data.positions - data.positions[0]) * move_factor
    xmin = float(data.x[0] + np.min(offset) - 0.12 * length)
    xmax = float(data.x[-1] + np.max(offset) + 0.12 * length)
    top = max(float(np.max(displayed_depths) * 1.20), float(np.max(data.depths) * 1.1))
    info_energy = data.energies - data.energies[0]
    expected_energy = data.input_work - data.viscous_loss
    error = data.balance_defect
    max_defect = float(np.max(np.abs(error)))

    colors = dict(bg="#071827", panel="#102b3d", grid="#365064", aqua="#16c8e8",
                  fill="#0b87b0", light="#e9f5fb", subtle="#93bed0", gold="#ffd179")
    fig = plt.figure(figsize=(12.8, 7.2), facecolor=colors["bg"])
    fig.text(.055, .946, "S L O S H I N G   T A N K", color=colors["light"],
             fontsize=22, fontweight="bold", family="DejaVu Sans")
    scheme_label = "Crank–Nicolson" if data.time_scheme == "cn" else "BDF2"
    fig.text(.057, .904,
             f"SCRIMP 1.2.0  |  1D viscous shallow water  |  {scheme_label}",
             color=colors["subtle"], fontsize=11)

    ax = fig.add_axes([.055, .39, .89, .47], facecolor=colors["panel"])
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-0.18 * top, 1.08 * top)
    ax.axis("off")
    ax.plot([xmin, xmax], [-0.09 * top, -0.09 * top],
            color=colors["grid"], lw=2.0, zorder=1)
    ax.text(xmin + .025 * length, 1.02 * top,
            "LABORATORY FRAME  ·  MOTION AMPLIFIED FOR DISPLAY",
            color=colors["subtle"], fontsize=9)
    tank_outline, = ax.plot([], [], lw=6, color=colors["light"],
                            solid_capstyle="round", zorder=7)
    water_patch = Polygon(np.zeros((len(data.x) + 2, 2)), closed=True,
                          facecolor=colors["fill"], alpha=.92, zorder=2)
    ax.add_patch(water_patch)
    water_line, = ax.plot([], [], color=colors["aqua"], lw=3, zorder=5)
    rest_line, = ax.plot([], [], color=colors["gold"], lw=1,
                         ls="--", alpha=.7, zorder=4)
    time_badge = fig.text(.73, .92, "", color=colors["gold"], fontsize=18, ha="left")
    position_badge = fig.text(.055, .350, "", color=colors["subtle"], fontsize=10)
    scale_badge = fig.text(.945, .350, "", color=colors["gold"], fontsize=10, ha="right")

    def setup_chart(box, title, ylabel, datasets):
        chart = fig.add_axes(box, facecolor=colors["panel"])
        stacked = []
        for values, label, line_color, linestyle in datasets:
            chart.plot(data.times, values, color=line_color, lw=2,
                       ls=linestyle, label=label)
            stacked.append(np.asarray(values))
        chart.set_xlim(data.times[0], data.times[-1])
        combined = np.concatenate(stacked)
        d = float(np.ptp(combined))
        if d < 1e-11:
            d = max(1e-8, float(np.max(np.abs(combined))) * .15)
        chart.set_ylim(np.min(combined) - .16*d, np.max(combined) + .16*d)
        chart.grid(alpha=.22, color=colors["subtle"])
        chart.tick_params(colors=colors["subtle"], labelsize=8)
        chart.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
        for spine in chart.spines.values():
            spine.set_color(colors["grid"])
        chart.set_title(title, color=colors["light"], fontsize=11,
                        loc="left", pad=10)
        chart.set_xlabel("time [s]", color=colors["subtle"], fontsize=9)
        chart.set_ylabel(ylabel, color=colors["subtle"], fontsize=9)
        if len(datasets) > 1:
            chart.legend(loc="best", fontsize=8, frameon=False,
                         labelcolor=colors["light"])
        cursor = chart.axvline(data.times[0], color=colors["gold"],
                               ls="--", lw=1.5)
        marker, = chart.plot([data.times[0]], [datasets[0][0][0]], "o",
                              color=colors["gold"], ms=5)
        return cursor, marker

    v_cursor, v_marker = setup_chart(
        [.055, .117, .254, .174], "Tank velocity", "v [m/s]",
        [(data.velocities, "v", colors["aqua"], "-")])
    e_cursor, e_marker = setup_chart(
        [.377, .117, .29, .174], "Power balance (integrated)", "energy [J]",
        [(info_energy, "H - H0", colors["gold"], "-"),
         (expected_energy, "Work - losses", colors["aqua"], "--")])
    r_cursor, r_marker = setup_chart(
        [.738, .117, .205, .174], "Balance defect", "residual [J]",
        [(error, "defect", "#fa9289", "-")])
    fig.text(
        .5, .025,
        f"POWER BALANCE: ΔH = W_ext − D_visc + ε_num    |    max |ε_num| = {max_defect:.2e} J",
        color=colors["subtle"], fontsize=10, ha="center")

    n_frames = max(2, int(round(seconds * fps)))
    frame_times = np.linspace(float(data.times[0]), float(data.times[-1]), n_frames)

    def draw(k):
        t = float(frame_times[k])
        depth, pos, velocity, energy = interpolate_frame(data, t)
        displayed_pos = (pos - data.positions[0]) * move_factor
        local_height = reference_height + surface_scale * (depth - reference_height)
        xs = data.x + displayed_pos
        poly = np.column_stack((np.r_[xs[0], xs, xs[-1]],
                                np.r_[0.0, local_height, 0.0]))
        water_patch.set_xy(poly)
        water_line.set_data(xs, local_height)
        rest_line.set_data([xs[0], xs[-1]], [reference_height] * 2)
        tank_outline.set_data([xs[0], xs[0], xs[-1], xs[-1]],
                              [top, 0, 0, top])
        time_badge.set_text(f"t = {t:6.3f} s")
        position_badge.set_text(
            f"physical tank position = {pos:+.4e} m     |     velocity = {velocity:+.4e} m/s")
        scale_badge.set_text(
            f"surface ×{surface_scale:g}  ·  tank motion ×{move_factor:.1f}")
        for cursor, marker, value in (
            (v_cursor, v_marker, velocity),
            (e_cursor, e_marker, energy-data.energies[0]),
            (r_cursor, r_marker, float(np.interp(t, data.times, error))),
        ):
            cursor.set_xdata([t, t])
            marker.set_data([t], [value])
        return (water_patch, water_line, rest_line, tank_outline,
                time_badge, position_badge, scale_badge,
                v_cursor, v_marker, e_cursor, e_marker,
                r_cursor, r_marker)

    writer = animation.FFMpegWriter(
        fps=fps, codec="libx264", bitrate=2400,
        extra_args=["-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                    "-pix_fmt", "yuv420p", "-movflags", "+faststart"])
    try:
        with writer.saving(fig, str(output), dpi=dpi):
            for k in range(n_frames):
                draw(k)
                writer.grab_frame()
    finally:
        plt.close(fig)
    return output


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path("outputs/sloshing_tank.mp4"))
    parser.add_argument("--mesh", type=float, default=0.05)
    parser.add_argument("--time-step", type=float, default=0.002)
    parser.add_argument("--final-time", type=float, default=0.2)
    parser.add_argument("--seconds", type=float, default=10.0)
    parser.add_argument("--fps", type=int, default=25)
    parser.add_argument("--samples", type=int, default=250)
    parser.add_argument("--surface-scale", type=float, default=8.0,
                        help="Visual gain applied to the surface deviation")
    parser.add_argument("--balance-csv", type=Path, default=None,
                        help="CSV sidecar path (default: output stem + .balance.csv)")
    parser.add_argument("--scheme", choices=["cn", "bdf"], default="cn",
                        help="Time integration: Crank-Nicolson by default, or BDF2")
    parser.add_argument("--motion-scale", default="auto",
                        help="'auto' or a positive scale (1 = physical movement)")
    args = parser.parse_args(argv)

    # Set the headless backend *before* SCRIMP imports pyplot.
    import matplotlib
    matplotlib.use("Agg")
    from sloshing_tank import sloshing_tank

    tank = sloshing_tank(mesh_size=args.mesh, final_time=args.final_time,
                        time_step=args.time_step, save_step=args.time_step,
                        adaptive=False, plot=False, time_scheme=args.scheme)
    data = extract_movie_data(tank, samples=args.samples)
    output = render_movie(data, args.output, fps=args.fps,
                          seconds=args.seconds, surface_scale=args.surface_scale,
                          motion_scale=args.motion_scale)
    csv_path = args.balance_csv if args.balance_csv else output.with_suffix(".balance.csv")
    csv_output = save_balance_csv(data, csv_path)
    print(f"MP4 created: {output.resolve()}")
    print(f"Balance time histories: {csv_output.resolve()}")
    print(f"Max absolute energy-balance defect: {np.max(np.abs(data.balance_defect)):.6e}")


if __name__ == "__main__":
    main()
