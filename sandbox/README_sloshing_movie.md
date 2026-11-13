# Sloshing-tank MP4 visualisation (post-1.2.0 prototype)

This tool renders an **actual SCRIMP simulation**, not a prescribed wave
animation. The source model remains the **experimental** 1D viscous
shallow-water / moving-tank model from SCRIMP v1.2.0.

Run from the repository root, in the working SCRIMP Conda environment:

```bash
conda install -c conda-forge ffmpeg  # only if not yet installed
python sandbox/animate_sloshing_tank.py \
    --mesh 0.05 --time-step 0.002 --final-time 15 \
    --seconds 15 --fps 25 --scheme cn \
    --output sandbox/outputs/sloshing_tank.mp4
```

The resulting H.264 MP4 is saved to the existing git-ignored `outputs/`
folder. Matplotlib's **Agg** backend is used: no Qt window or Wayland/X11
display server is required. Check that `ffmpeg -version` works.

## What is shown

- The tank moves horizontally according to the global state `pos(t)`;
  its **physical** position and velocity appear in numeric readouts.
- The free surface is the computed finite-element water depth
  **`h(x,t)=a_h(x,t)/b`**, sampled using GetFEM interpolation, not an
  assumed sinusoidal profile.
- Three synchronized panels display tank velocity, **total energy change
  versus integrated input work minus viscous dissipation**, and the
  **energy/power-balance defect** `ε(t)`.
- Intermediate video frames are **linearly interpolated** between
  numerical solution snapshots; they are not new solver time steps.

For visibility, the default animation magnifies deviations of the water
surface from its **initial mean by a factor of 8**, and chooses an
automatic magnification of the tank displacement. The factors appear
**on screen**. Numeric values and diagnostic plots are never scaled.

To remove all visual exaggeration, use
`--surface-scale 1 --motion-scale 1`.

## Options and testing

```bash
python sandbox/animate_sloshing_tank.py --help
python -m unittest -v tests.test_sloshing_movie
```

- `--mesh`, `--time-step`, `--final-time`: actual simulation settings.
- `--seconds`, `--fps`: MP4 playback settings, independently of simulation time.
- `--samples`: number of spatial interpolation points.
- `--surface-scale`: display gain applied only to the perturbation of the depth.
- `--motion-scale`: `auto` or positive numeric value; `1` means physical scale.
- `--output`: MP4 destination.

The seven renderer/CSV/FFmpeg tests passed without a GetFEM installation.
**A complete integration test with a real SCRIMP/GetFEM simulation remains
necessary on the numerical workstation.**

**Scientific caution:** a successful animation does not validate the
tank's impermeable-wall conditions, energy accounting or spatial
convergence. The model remains an experimental sandbox prototype.

## Crank–Nicolson and the power balance

The simulation now uses **Crank–Nicolson** (`--scheme cn`) by default.
To reproduce earlier BDF2 experiments explicitly, use `--scheme bdf`.
The time step is held fixed for movie generation.

In the experimental mathematical model, with tank velocity `v`,
external force `F(t)`, fluid cross-section `a_h`, and fluid relative
co-energy `e_m`, the nominal continuous power law (assuming closed,
impermeable boundaries) is

```text
dH/dt = F(t) * v(t) - D_visc(t)
D_visc = integral(rho * nu * a_h * |d_x(e_m)|^2) dx  >= 0.
```

From the **saved numerical solution**, the script computes cumulative
work and dissipation using the **trapezoidal rule**:

```text
W(t) = integral_0^t F(s)*v(s) ds
D(t) = integral_0^t D_visc(s) ds
epsilon(t) = H(t) - H(0) - [W(t) - D(t)].
```

The MP4 compares `H(t)-H(0)` to `W(t)-D(t)`, and plots
`epsilon(t)` separately with its maximum absolute value. A near-zero
defect is a **measured** numerical property, not an assumption.
Crank–Nicolson is second-order and time-symmetric, but for this nonlinear
system it is **not in general an exactly symplectic or
energy-preserving discretization**; the postprocessed trapezoidal power
integral also introduces quadrature error.

Each run also writes a **CSV sidecar**, e.g.
`sandbox/outputs/sloshing_tank.balance.csv`, with the original saved
timestamps, Hamiltonian, power input, viscous dissipation, cumulative
work and losses, expected energy increment, and balance defect.
All these quantities are **unscaled physical values**. The movie's
surface and displacement amplification applies only to its depiction.

A successful small energy defect does **not** validate the omitted wall
constraints or spatial convergence of the experimental PDE model.
