# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Headless checks for the experimental sloshing-tank movie renderer."""

import tempfile
from pathlib import Path
import unittest

import numpy as np

from sandbox.animate_sloshing_tank import (
    MovieData, interpolate_frame, render_movie, validate_data,
)


class TestSloshingMovie(unittest.TestCase):
    def setUp(self):
        self.data = MovieData(
            times=np.array([0.0, 0.5, 1.0]),
            x=np.linspace(0.0, 2.0, 9),
            depths=np.array([np.ones(9) * 2.0,
                             np.ones(9) * 2.1,
                             np.ones(9) * 2.0]),
            positions=np.array([0., 0.01, 0.02]),
            velocities=np.array([0., .02, .0]),
            energies=np.array([40., 39.9, 39.8]),
        )

    def test_saved_state_interpolation(self):
        depth, pos, vel, ham = interpolate_frame(self.data, 0.25)
        np.testing.assert_allclose(depth, 2.05)
        self.assertAlmostEqual(pos, 0.005)
        self.assertAlmostEqual(vel, 0.01)
        self.assertAlmostEqual(ham, 39.95)

    def test_negative_depth_is_rejected(self):
        self.data.depths[1, 3] = -0.5
        with self.assertRaisesRegex(ValueError, "non-positive"):
            validate_data(self.data)

    def test_duplicate_times_are_rejected(self):
        self.data.times[1] = 0.0
        with self.assertRaisesRegex(ValueError, "increasing"):
            validate_data(self.data)

    def test_mp4_rendering(self):
        from matplotlib import animation
        if not animation.writers.is_available("ffmpeg"):
            self.skipTest("ffmpeg is unavailable")
        with tempfile.TemporaryDirectory() as tmp:
            output = render_movie(self.data, Path(tmp) / "movie.mp4",
                                  fps=3, seconds=1.0, dpi=48)
            self.assertTrue(output.is_file())
            self.assertGreater(output.stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
