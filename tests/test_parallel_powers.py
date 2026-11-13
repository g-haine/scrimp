# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Regression tests for rank-safe Hamiltonian/power plotting (issue #6)."""

import unittest
from unittest.mock import MagicMock, call, patch

import numpy as np

import scrimp.dphs as dphs_module
from scrimp.dphs import DPHS
from scrimp.hamiltonian import Hamiltonian
from scrimp.port import Port
import scrimp.port as port_module


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
                with patch.object(dphs_module, "rank", rank), \
                     patch.object(dphs_module.plt, "figure", return_value=fig) as figure, \
                     patch.object(dphs_module.plt, "show") as show, \
                     patch.object(DPHS, "plot_powers", autospec=True) as plot_powers:
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
                with patch.object(dphs_module, "rank", rank), \
                     patch.object(dphs_module.plt, "figure", return_value=MagicMock()), \
                     patch.object(dphs_module.plt, "show"), \
                     patch.object(DPHS, "plot_powers", autospec=True) as plot_powers:
                    self.system.plot_Hamiltonian(with_powers=False)
                    plot_powers.assert_not_called()


class TestPortPowerAssembly(unittest.TestCase):
    def test_getfem_time_is_synchronized_for_each_power_sample(self):
        port = Port("test", "f", "e", "scalar-field",
                    algebraic=True, substituted=False)
        model = MagicMock()
        domain = MagicMock()
        domain._int_method = {0: "integration"}
        solution = {
            "t": [0.0, 0.1, 0.2],
            "z": [
                np.array([0.0, 2.0]),
                np.array([2.0, 4.0]),
                np.array([4.0, 6.0]),
            ],
        }
        with patch.object(port_module.gf, "asm", side_effect=[1.0, 2.0, 3.0]):
            port.compute(solution, model, domain)

        self.assertEqual(model.set_time.call_args_list,
                         [call(0.0), call(0.1), call(0.2)])
        np.testing.assert_array_equal(
            model.to_variables.call_args_list[1].args[0], [1.0, 3.0]
        )
        np.testing.assert_array_equal(
            model.to_variables.call_args_list[2].args[0], [3.0, 5.0]
        )
        self.assertEqual(port.get_power(), [1.0, 2.0, 3.0])


if __name__ == "__main__":
    unittest.main()
