# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Tests for Hamiltonian and Term as implemented in SCRIMP."""

import unittest

from scrimp.hamiltonian import Hamiltonian, Term


class TestHamiltonian(unittest.TestCase):
    def setUp(self):
        self.hamiltonian = Hamiltonian("Total mechanical energy")
        self.term = Term("Kinetic energy", "p*p/(2*m)", [1])

    def test_name_and_rename(self):
        self.assertEqual(self.hamiltonian.get_name(), "Total mechanical energy")
        self.hamiltonian.set_name("New energy")
        self.assertEqual(self.hamiltonian.get_name(), "New energy")

    def test_initially_empty(self):
        self.assertEqual(len(self.hamiltonian), 0)
        self.assertEqual(self.hamiltonian.get_terms(), [])
        self.assertFalse(self.hamiltonian.get_is_computed())

    def test_register_hamiltonian_terms(self):
        self.hamiltonian.add_term(self.term)
        self.assertEqual(len(self.hamiltonian), 1)
        self.assertIs(self.hamiltonian[0], self.term)
        self.assertEqual(list(self.hamiltonian), [self.term])

    def test_get_terms_returns_copy(self):
        self.hamiltonian.add_term(self.term)
        terms = self.hamiltonian.get_terms()
        terms.clear()
        self.assertEqual(self.hamiltonian.get_terms(), [self.term])

    def test_add_term_rejects_unsupported_object(self):
        with self.assertRaises(AssertionError):
            self.hamiltonian.add_term("not a Term")

    def test_computation_flag(self):
        self.assertFalse(self.hamiltonian.get_is_computed())
        self.hamiltonian.set_is_computed()
        self.assertTrue(self.hamiltonian.get_is_computed())


if __name__ == "__main__":
    unittest.main()
