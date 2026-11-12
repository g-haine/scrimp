# Sloshing-tank MP4 visualisation (post-1.2.0 prototype)

This tool renders an **actual SCRIMP simulation**, not a prescribed wave
animation. The source model remains the **experimental** 1D viscous
shallow-water / moving-tank model from SCRIMP v1.2.0.

Run from the repository root, in the working SCRIMP Conda environment:

```bash
conda install -c conda-forge ffmpeg  # only if not yet installed
python sandbox/animate_sloshing_tank.py \
    --mesh 0.05 --time-step 0.002 --final-time 0.2 \
    --seconds 10 --fps 25 --output outputs/sloshing_tank.mp4
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
- Two synchronized plots display the physical tank velocity and the
  **change in total Hamiltonian** `H(t)-H(0)`.
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

The pure renderer and FFmpeg smoke tests passed without a GetFEM installation.
**A complete integration test with a real SCRIMP/GetFEM simulation remains
necessary on the numerical workstation.**

**Scientific caution:** a successful animation does not validate the
tank's impermeable-wall conditions, energy accounting or spatial
convergence. The model remains an experimental sandbox prototype.
