# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Regression tests for rank-safe Hamiltonian/power plotting (issue #6)."""

import unittest
from unittest.mock import MagicMock, patch

import numpy as np

import scrimp.dphs as dphs_module
from scrimp.dphs import DPHS
from scrimp.hamiltonian import Hamiltonian


class TestParallelPowerPlot(unittest.TestCase):
    def setUp(self):
        # No GetFEM or PETSc solver is needed to test MPI call symmetry.
        self.system = object.__new__(DPHS)
        self.system.solution = {"t": [0.0, 1.0]}
        self.system.hamiltonian = Hamiltonian("Energy")
        self.system.hamiltonian.set_is_computed()

    def test_plot_powers_is_entered_on_every_rank(self):
        for rank in (0, 1):
            with self.subTest(rank=rank):
                fig = MagicMock()
                ax = fig.add_subplot.return_value
                with (
                    patch.object(dphs_module, "rank", rank),
                    patch.object(dphs_module.plt, "figure", return_value=fig) as figure,
                    patch.object(dphs_module.plt, "show") as show,
                    patch.object(DPHS, "plot_powers", autospec=True) as plot_powers,
                ):
                    self.system.plot_Hamiltonian(with_powers=True)
                    plot_powers.assert_called_once()
                    args, kwargs = plot_powers.call_args
                    self.assertIs(args[0], self.system)
                    self.assertIs(args[1], ax if rank == 0 else None)
                    np.testing.assert_array_equal(kwargs["HamTot"], [0.0, 0.0])
                    self.assertEqual(figure.call_count, 1 if rank == 0 else 0)
                    self.assertEqual(show.call_count, 1 if rank == 0 else 0)

    def test_without_powers_only_root_renders(self):
        for rank in (0, 1):
            with self.subTest(rank=rank):
                with (
                    patch.object(dphs_module, "rank", rank),
                    patch.object(dphs_module.plt, "figure", return_value=MagicMock()),
                    patch.object(dphs_module.plt, "show"),
                    patch.object(DPHS, "plot_powers", autospec=True) as plot_powers,
                ):
                    self.system.plot_Hamiltonian(with_powers=False)
                    plot_powers.assert_not_called()


if __name__ == "__main__":
    unittest.main()
