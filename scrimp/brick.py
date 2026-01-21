# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             brick.py
- authors:          Giuseppe Ferraro, Ghislain Haine
- date:             31 may 2023
- brief:            class for brick object
"""

from typing import Any, Dict


class Brick:
    """This class defines a Brick."""

    def __init__(
        self,
        name: str,
        form: str,
        regions: list,
        linear: bool = True,
        dt: bool = False,
        position: str = "constitutive",
        explicit: bool = False,
        mesh_id: int = 0,
        structure: Any | None = None,
    ):
        """_summary_

        Args:
            name (str): the name of the brick, will be used mainly for plotting purpose.
            form (str): the form in GWFL getfem language.
            regions (list): the regions of mesh where the form applies.
            linear (bool, optional): parameter to help easy identification of linear bricks. Defaults to True.
            dt (bool, optional): parameter to help easy identification of matrices applied to time-derivative of a variable (e.g. mass matrices). Defaults to False.
            position (str, optional): parameter to help easy identification of "where" is the form: in the Dirac structure ('flow' side or 'effort' side), or in the 'constitutive' relations. This serves for both the time-resolution and plotting purposes. Defaults to "constitutive".
            explicit (str, optional): parameter to indicate wether the brick has to be treated explicitly in the time-stepping scheme.
            mesh_id (int, optional): the id of the mesh where the form applies. Defaults to 0.
        """

        self._name = name
        self._id_bricks = []
        self._form = form
        self._mesh_id = mesh_id
        self._regions = regions
        self._linear = linear
        self._dt = dt
        self._position = position
        self._explicit = explicit
        self._structure = structure

    def add_id_brick_to_list(self, id_brick: int):
        """This function adds a brick ID to the brick ID list.

        Args:
            id_brick (int): the id of the brick
        """

        self._id_bricks.append(id_brick)

    def get_name(self) -> str:
        """This function returns the name of the brick.

        Returns:
            str: the name of the brick, will be used mainly for plotting purpose
        """

        return self._name

    def get_id_bricks(self) -> list:
        """This function returns the list of integers related to the ids of the bricks.

        Returns:
            list: the list of integers related to the ids of the bricks.
        """

        return self._id_bricks.copy()

    def get_form(self) -> str:
        """This function returns the form of the brick.

        Returns:
            str: the form in GWFL getfem language.
        """

        return self._form

    def get_mesh_id(self) -> int:
        """This function returns the ID of the brick.

        Returns:
            int: the id of the mesh where the form applies.
        """

        return self._mesh_id

    def get_regions(self) -> list:
        """This function returns the regions of the brick.

        Returns:
            list: the regions of mesh where the form applies.
        """

        return self._regions

    def get_linear(self) -> bool:
        """This function returns the boolean that defines wether the brick is linear or not.

        Returns:
            bool: parameter to help easy identification of linear bricks.
        """

        return self._linear

    def get_dt(self) -> bool:
        """This function returns the boolean that defines wether the matrices applied to time-derivative of a variable or not.

        Returns:
            bool: parameter to help easy identification of matrices applied to time-derivative of a variable (e.g. mass matrices).
        """

        return self._dt

    def get_position(self) -> int:
        """This function returns the id of the position where the form of the brick applies.

        Returns:
            int: the id of the mesh where the form applies. Defaults to 0.
        """

        return self._position

    def get_explicit(self) -> bool:
        """This function returns the boolean that defines wether the brick is explicit or not.

        Returns:
            bool: parameter to help easy identification of explicit bricks.
        """

        return self._explicit

    def disable_id_bricks(self, gf_model):
        """This function disable the brick in the getfem model."""
        gf_model.disable_bricks(self._id_bricks)

    def enable_id_bricks(self, gf_model):
        """This function enable the brick in the getfem model."""
        gf_model.enable_bricks(self._id_bricks)

    def get_structure(self) -> Any | None:
        """Return the symbolic structure associated with the brick."""

        return self._structure

    def attach_structure(self, structure: Any) -> None:
        """Attach a symbolic constitutive relation to the brick."""

        self._structure = structure

    def substitutions(self) -> Dict[str, Any]:
        """Return constitutive substitutions derived from the structure."""

        from scrimp.structure.core import ConstitutiveRelation

        if isinstance(self._structure, ConstitutiveRelation):
            mapping = self._structure.substitution_map()
            return {str(symbol): expr for symbol, expr in mapping.items()}
        return {}
