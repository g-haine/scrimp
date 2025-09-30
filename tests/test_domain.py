import unittest

from scrimp import Domain
from scrimp.io.schema_loader import DomainSchema, IntegrationRuleSchema, MeshSchema


class TestDomain(unittest.TestCase):
    def test_legacy_initialisation(self):
        domain = Domain("Interval", {"L": 2.0, "h": 0.15})
        self.assertEqual(domain.get_name(), "Interval")
        self.assertTrue(domain.get_isSet())
        self.assertIn("Interval", domain.get_mesh_labels())

    def test_schema_initialisation(self):
        schema = DomainSchema(
            name="IntervalDomain",
            meshes=[
                MeshSchema(
                    id="omega",
                    source="builtin",
                    generator="Interval",
                    parameters={"L": 2.0, "h": 0.1},
                    integration=IntegrationRuleSchema(family="gauss", order=3),
                )
            ],
        )
        domain = Domain(schema)
        self.assertEqual(domain.get_name(), "IntervalDomain")
        self.assertListEqual(domain.get_mesh_labels(), ["omega"])
        self.assertEqual(domain.get_dim()[0], 1)
        exported = domain.get_schema()
        self.assertIsNotNone(exported)
        self.assertEqual(exported.dict(), schema.dict())


if __name__ == "__main__":
    unittest.main()
