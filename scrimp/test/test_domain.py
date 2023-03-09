import unittest
from src.domain import Domain
import utils.mesh


class TestDomain(unittest.TestCase):
    def test_name(self):
        domain = Domain("Rectangle", {"L": 2.0, "l": 1.0, "h": 0.15})
        self.assertEqual(domain.get_name(), "Rectangle")  # add assertion here

    def test_builtins_methods(self):
        built_in_methods = dir(utils.mesh)

        name = "Rectangle"
        # print(built_in_methods)
        self.assertTrue(name in built_in_methods)

        name = "ciao"
        self.assertFalse(name in built_in_methods)

    def test_mesh(self):
        name = "Rectangle"
        parameters = {"L": 2.0, "l": 1.0, "h": 0.15}
        domain = Domain(name, parameters)
        built_in_methods = dir(utils)
        gf_mesh = None

        if name in built_in_methods:
            gf_mesh = utils.mesh.Rectangle(parameters, 0)
            self.assertListEqual(domain.get_mesh(), gf_mesh[0])  # add assertion here

        else:
            self.assertIsNone(gf_mesh)

    def test_dim(self):
        name = "Rectangle"
        parameters = {"L": 2.0, "l": 1.0, "h": 0.15}
        domain = Domain(name, parameters)
        built_in_methods = dir(utils)
        gf_mesh = None

        if name in built_in_methods:
            gf_mesh = utils.mesh.Rectangle(parameters, 0)
            self.assertListEqual(domain.get_dim(), gf_mesh[1])  # add assertion here

        else:
            self.assertIsNone(gf_mesh)

    def test_subdomains(self):
        name = "Rectangle"
        parameters = {"L": 2.0, "l": 1.0, "h": 0.15}
        domain = Domain(name, parameters)
        built_in_methods = dir(utils)
        gf_mesh = None

        if name in built_in_methods:
            gf_mesh = utils.mesh.Rectangle(parameters, 0)
            self.assertListEqual(
                domain.get_subdomains(), gf_mesh[2]
            )  # add assertion here

        else:
            self.assertIsNone(gf_mesh)

    def test_boundaries(self):
        name = "Rectangle"
        parameters = {"L": 2.0, "l": 1.0, "h": 0.15}
        domain = Domain(name, parameters)
        built_in_methods = dir(utils)
        gf_mesh = None

        if name in built_in_methods:
            gf_mesh = utils.mesh.Rectangle(parameters, 0)
            self.assertListEqual(
                domain.get_boundaries(), gf_mesh[3]
            )  # add assertion here

        else:
            self.assertIsNone(gf_mesh)

    def test_isSet(self):
        name = "Rectangle"
        parameters = {"L": 2.0, "l": 1.0, "h": 0.15}
        domain = Domain(name, parameters)
        self.assertTrue(domain.get_isSet())

        name = "Ciao"
        parameters = {"L": 2.0, "l": 1.0, "h": 0.15}
        domain = Domain(name, parameters)
        self.assertFalse(domain.get_isSet())


if __name__ == "__main__":
    unittest.main()
