import unittest
from scrimp.port import Port, Parameter
from scrimp.structure import (
    LagrangeSubspace,
    DiracStructure,
    ConstitutiveRelation,
)


class TestPort(unittest.TestCase):
    def test_name(self):
        port = Port("name", "flow", "effort", "kind", 0, True, True, 1)
        self.assertEqual(port.get_name(), "name")  # add assertion here

    def test_flow(self):
        port = Port("name", "flow", "effort", "kind", 0, True, True, 1)
        self.assertEqual(port.get_flow(), "flow")  # add assertion here

    def test_effort(self):
        port = Port("name", "flow", "effort", "kind", 0, True, True, 1)
        self.assertEqual(port.get_effort(), "effort")  # add assertion here

    def test_kind(self):
        port = Port("name", "flow", "effort", "kind", 0, True, True, 1)
        self.assertEqual(port.get_kind(), "kind")  # add assertion here

    def test_mesh_id(self):
        port = Port("name", "flow", "effort", "kind", 0, True, True, 1)
        self.assertEqual(port.get_mesh_id(), 0)  # add assertion here

    def test_algebraic(self):
        port = Port("name", "flow", "effort", "kind", 0, True, True, 1)
        self.assertEqual(port.get_algebraic(), True)  # add assertion here

    def test_substituted(self):
        port = Port("name", "flow", "effort", "kind", 0, True, True, 1)
        self.assertEqual(port.get_substituted(), True)  # add assertion here

    def test_region(self):
        port = Port("name", "flow", "effort", "kind", 0, True, True, 1)
        self.assertEqual(port.get_region(), 1)  # add assertion here

    def test_set_fem(self):
        pass
        # port = Port("name", "flow", "effort", "kind", 0, True, True, 1)
        # # try to set a coport when the port doesn't have one yet
        # coport1 = Coport("name_coport1", "descrtiption_coport", port, True)
        # port.set_coport(coport1)
        # self.assertEqual(coport1, port.get_coport())
        #
        # # try to set a coport when the port does have one already
        # coport2 = Coport("name_coport2", "descrtiption_coport", port, True)
        # port.set_coport(coport2)
        # self.assertEqual(coport1, port.get_coport())

    def test_fem(self):
        port = Port("name", "flow", "effort", "kind", 0, True, True, 1)
        self.assertEqual(port.get_fem(), None)  # add assertion here

    def test_get_parameter(self):
        port = Port("name", "flow", "effort", "kind", 0, True, True, 1)
        self.assertIsNone(port.get_parameter("not_present_name"))

        parameter = Parameter(
            "name_para", "description", "kind_param", "expression", port.get_name()
        )
        port.add_parameter(parameter)
        self.assertEqual(port.get_parameter("name_para"), parameter)

    def test_add_parameter(self):
        port = Port("name", "flow", "effort", "kind", 0, True, True, 1)
        self.assertFalse(port.add_parameter(None))

        parameter = Parameter(
            "name_para", "description", "kind_param", "expression", port.get_name()
        )
        self.assertTrue(port.add_parameter(parameter))
        self.assertEqual(port.get_parameter("name_para"), parameter)

    # TODO add fem obgect
    # def test_init_parameter(self):
    #     port = Port("name", "flow", "effort", "kind", 0, True, True, 1)
    #     parameter = Parameter(
    #         "name_para", "description", "kind_param", "expression", port.get_name()
    #     )
    #     port.add_parameter(parameter)
    #     self.assertEqual(port, port.init_parameter("name_para","expression"))
    #

    def test_str(self):
        port = Port("name", "flow", "effort", "kind", 0, True, True, 1)
        out = port.__str__()
        print(out)
        self.assertEqual(
            out,
            f"{port.get_name()}, {port.get_flow()}, {port.get_effort()}, {port.get_kind()}, {str(port.get_mesh_id())}, {str(port.get_algebraic())}, {port.get_parameters()}",
        )

    def test_structure_validation(self):
        subspace = LagrangeSubspace("f", "e", metadata={"regions": [1]})
        port = Port("name", "flow", "effort", "kind", structure=subspace)
        cancelling = LagrangeSubspace(
            "f", "e", metadata={"regions": [1], "sign": -1}
        )
        dirac = DiracStructure([subspace, cancelling])

        port.validate_power_balance(dirac)

    def test_constitutive_substitutions(self):
        relation = ConstitutiveRelation("Ohm", {"e": "R*f"})
        port = Port("name", "flow", "effort", "kind", structure=relation)

        subs = port.substitutions()
        self.assertIn("e", subs)


if __name__ == "__main__":
    unittest.main()
