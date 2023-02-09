import unittest
from scrimp.dphs.port import Port, Parameter


class TestParameter(unittest.TestCase):

    def test_name(self):
        parameter = Parameter("name", "description", "kind", "expression", "name_port")
        self.assertEqual(parameter.get_name(), "name")  # add assertion here

    def test_description(self):
        parameter = Parameter("name", "description", "kind", "expression", "name_port")
        self.assertEqual(parameter.get_description(), "description")  # add assertion here

    def test_kind(self):
        parameter = Parameter("name", "description", "kind", "expression", "name_port")
        self.assertEqual(parameter.get_kind(), "kind")  # add assertion here

    def test_expression(self):
        parameter = Parameter("name", "description", "kind", "expression", "name_port")
        self.assertEqual(parameter.get_expression(), "expression")  # add assertion here

    def test_name_port(self):
        parameter = Parameter("name", "description", "kind", "expression", "name_port")
        self.assertEqual(parameter.get_name_port(), "name_port")  # add assertion here


if __name__ == '__main__':
    unittest.main()
