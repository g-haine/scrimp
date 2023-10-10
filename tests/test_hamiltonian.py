import unittest
from scrimp import Term, Hamiltonian

class TestHamiltonian(unittest.TestCase):
    def test_add_term(self):
        hamiltonian = Hamiltonian("hamiltonian")
        term = Term("Potential energy", "0.5*q.T.q", [1])
        hamiltonian.add_term(term)

        self.assertEqual(hamiltonian.get_terms()[0], term)

    def test_get_name(self):
        hamiltonian = Hamiltonian("hamiltonian")
        self.assertEqual(hamiltonian.get_name(), "hamiltonian")

    def test_set_name(self):
        hamiltonian = Hamiltonian("hamiltonian")
        name = "test_name"
        hamiltonian.set_name(name)
        self.assertEqual(hamiltonian.get_name(), name)

    def test_get_terms(self):
        hamiltonian = Hamiltonian("hamiltonian")
        self.assertListEqual(hamiltonian.get_terms(), [])

    def test_set_is_computed(self):
        hamiltonian = Hamiltonian("hamiltonian")
        hamiltonian.set_is_computed()
        self.assertTrue(hamiltonian.get_is_computed())

    def test_get_is_computed(self):
        hamiltonian = Hamiltonian("hamiltonian")
        self.assertFalse(hamiltonian.get_is_computed())
        hamiltonian.set_is_computed()
        self.assertTrue(hamiltonian.get_is_computed())

if __name__ == "__main__":
    unittest.main()
