import unittest
from scrimp.brick import Brick
from scrimp.structure import ConstitutiveRelation, ensure_brick


class TestBrick(unittest.TestCase):
    def test_name(self):
        brick = Brick("name", "form", [1, 2], True, False, "constitutive", 0)
        self.assertEqual(brick.get_name(), "name")

    def test_form(self):
        brick = Brick("name", "form", [1, 2], True, False, "constitutive", 0)
        self.assertEqual(brick.get_form(), "form")

    def test_id_bricks(self):
        brick = Brick("name", "form", [1, 2], True, False, "constitutive", 0)
        self.assertListEqual(brick.get_id_bricks(), [])

    def test_linear(self):
        brick = Brick("name", "form", [1, 2], True, False, "constitutive", 0)
        self.assertTrue(brick.get_linear())

    def test_dt(self):
        brick = Brick("name", "form", [1, 2], True, False, "constitutive", 0)
        self.assertFalse(brick.get_dt())

    def test_position(self):
        brick = Brick("name", "form", [1, 2], True, False, "constitutive", 0)
        self.assertEqual(brick.get_position(), "constitutive")

    def test_get_mesh_id(self):
        brick = Brick("name", "form", [1, 2], True, False, "constitutive", 0)
        self.assertEqual(brick.get_mesh_id(), 0)

    def test_symbolic_structure(self):
        relation = ConstitutiveRelation("Ohm", {"e": "R*f"})
        brick = Brick(
            "constitutive",
            "form",
            [1],
            structure=relation,
        )
        self.assertIs(brick.get_structure(), relation)
        subs = brick.substitutions()
        self.assertIn("e", subs)

    def test_ensure_brick(self):
        relation = ConstitutiveRelation(
            "Ohm", {"e": "R*f"}, metadata={"regions": [1], "mesh_id": 0}
        )
        brick = ensure_brick(relation)
        self.assertIs(brick.get_structure(), relation)


if __name__ == "__main__":
    unittest.main()
