import unittest
from scrimp import Domain, Term
import getfem as gf


class TestTerm(unittest.TestCase):
    def test_get_description(self):
        term = Term("Potential energy", "0.5*q.T.q", [1])
        self.assertEqual(term.get_description(), "Potential energy")

    def test_get_expression(self):
        term = Term("Potential energy", "0.5*q.T.q", [1])
        self.assertEqual(term.get_expression(), "0.5*q.T.q")

    def test_get_regions(self):
        term = Term("Potential energy", "0.5*q.T.q", [1])
        self.assertEqual(term.get_regions(), [1])

    def test_get_mesh_id(self):
        term = Term("Potential energy", "0.5*q.T.q", [1])
        self.assertEqual(term.get_mesh_id(), 0)

    def test_get_values(self):
        term = Term("Potential energy", "0.5*q.T.q", [1])
        self.assertEqual(term.get_values(), [])

    def test_set_value(self):
        domain = Domain("Rectangle", {"L": 2.0, "l": 1.0, "h": 0.15})
        gf_model = gf.Model("real")
        term = Term("Potential energy", "0.5*q.T.q", [1])
        term_value_at_t = 0.0
        term.set_value(term_value_at_t)
        self.assertListEqual(term.get_values(), [0.0])

if __name__ == "__main__":
    unittest.main()
