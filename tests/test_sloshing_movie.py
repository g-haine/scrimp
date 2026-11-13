# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Solver-independent checks for the sloshing-tank MP4 and energy diagnostics."""

import csv
import tempfile
from pathlib import Path
import unittest

import numpy as np

from sandbox.animate_sloshing_tank import (
    MovieData, interpolate_frame, render_movie, save_balance_csv,
    save_classic_hamiltonian, validate_data,
)


class TestSloshingMovie(unittest.TestCase):
    def setUp(self):
        self.data = MovieData(
            times=np.array([0.0, 0.5, 1.0]),
            x=np.linspace(0., 2., 9),
            depths=np.array([np.ones(9)*2., np.ones(9)*2.1, np.ones(9)*2.]),
            positions=np.array([0., .01, .02]),
            velocities=np.array([0., .02, .0]),
            energies=np.array([40., 39.9, 39.8]),
            input_work=np.zeros(3),
            viscous_loss=np.array([0., .1, .2]),
            balance_defect=np.zeros(3),
            input_power=np.zeros(3),
            viscous_power=np.ones(3)*.2,
            time_scheme="cn",
        )

    def test_interpolation(self):
        depth, pos, vel, energy = interpolate_frame(self.data, .25)
        np.testing.assert_allclose(depth, 2.05)
        self.assertAlmostEqual(pos, .005)
        self.assertAlmostEqual(vel, .01)
        self.assertAlmostEqual(energy, 39.95)

    def test_balance_identity(self):
        validate_data(self.data)
        np.testing.assert_allclose(
            self.data.energies - self.data.energies[0],
            self.data.input_work - self.data.viscous_loss
            + self.data.balance_defect,
        )

    def test_inconsistent_balance_is_rejected(self):
        self.data.balance_defect[-1] = .005
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            validate_data(self.data)

    def test_negative_depth_is_rejected(self):
        self.data.depths[1, 3] = -.5
        with self.assertRaisesRegex(ValueError, "non-positive"):
            validate_data(self.data)

    def test_duplicate_times_are_rejected(self):
        self.data.times[1] = .0
        with self.assertRaisesRegex(ValueError, "increasing"):
            validate_data(self.data)

    def test_csv_export_has_real_values(self):
        with tempfile.TemporaryDirectory() as temp:
            output = save_balance_csv(self.data, Path(temp) / "balance.csv")
            with output.open(encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 3)
            self.assertEqual(rows[0]["time_scheme"], "cn")
            self.assertAlmostEqual(float(rows[-1]["balance_defect"]), 0.)

    def test_native_scrimp_hamiltonian_image(self):
        """The native SCRIMP call is preserved, not redrawn from sampled data."""
        from matplotlib import pyplot as plt

        class FakeDPHS:
            def __init__(self, scratch):
                self.scratch = scratch
                self.calls = []

            def plot_Hamiltonian(self, with_powers=True, save_figure=False,
                                 filename="Hamiltonian.png"):
                self.calls.append((with_powers, save_figure, filename))
                fig, ax = plt.subplots()
                ax.plot([0, 1], [1, .9], label="Total Hamiltonian")
                ax.legend()
                fig.savefig(self.scratch / filename, dpi=70)
                # Real SCRIMP also calls plt.show(): Agg makes that harmless.
                plt.show()

        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            scratch = folder / "scratch" / "png"
            scratch.mkdir(parents=True)
            fake = FakeDPHS(scratch)
            for suffix, signature in ((".png", b"\\x89PNG"),
                                      (".jpg", b"\\xff\\xd8")):
                destination = folder / "video" / ("tank_hamiltonian" + suffix)
                returned = save_classic_hamiltonian(
                    fake, destination, source_dir=scratch
                )
                self.assertEqual(returned, destination)
                self.assertTrue(destination.is_file())
                self.assertTrue(destination.read_bytes().startswith(signature))
            self.assertEqual(len(fake.calls), 2)
            for with_powers, save_figure, filename in fake.calls:
                self.assertFalse(with_powers)
                self.assertTrue(save_figure)
                self.assertTrue(filename.startswith("_slosh_hamiltonian_"))
            self.assertFalse(list(scratch.iterdir()))

    def test_mp4_rendering(self):
        from matplotlib import animation
        if not animation.writers.is_available("ffmpeg"):
            self.skipTest("ffmpeg is unavailable")
        with tempfile.TemporaryDirectory() as temp:
            output = render_movie(
                self.data, Path(temp) / "movie.mp4",
                fps=3, seconds=1.0, dpi=48,
            )
            self.assertTrue(output.is_file())
            self.assertGreater(output.stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
