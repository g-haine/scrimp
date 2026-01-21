# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             domain.py
- authors:          Giuseppe Ferraro, Ghislain Haine
- date:             31 may 2023
- brief:            class for domain object
"""

import petsc4py
import sys
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

petsc4py.init(sys.argv)
from petsc4py import PETSc

comm = PETSc.COMM_WORLD
rank = comm.getRank()

import scrimp.utils.mesh
import getfem as gf
import logging

from scrimp.io.schema_loader import (
    DomainSchema,
    IntegrationRuleSchema,
    MeshSchema,
    load_domain_schema,
)

import scrimp.utils.config
outputs_path = scrimp.utils.config.outputs_path


class Domain:
    """A class handling meshes and indices for regions

    Lists are used to handle interconnection of pHs, allowing for several meshes in the dpHs.
    """

    def __init__(
        self,
        name: Union[str, DomainSchema],
        parameters: Optional[dict] = None,
        refine: int = 0,
        terminal: int = 1,
        schema: Optional[DomainSchema] = None,
    ):
        """Constructor of the `domain` member of a dpHs"""

        self._name = name
        self._isSet = False  #: A boolean to check if the domain has been set
        self._mesh: List[gf.Mesh] = []
        self._subdomains: List[Dict[str, int]] = []
        self._boundaries: List[Dict[str, int]] = []
        self._dim: List[int] = []
        self._int_method: List[gf.MeshIm] = []
        self._mesh_labels: List[str] = []
        self._mesh_schemas: List[Optional[MeshSchema]] = []

        parameters = parameters or {}

        if isinstance(name, DomainSchema):
            schema = name
        elif schema is None and isinstance(name, (str, Path)) and os.path.isfile(name):
            # Allow passing a schema file path directly
            try:
                schema = load_domain_schema(name)
            except Exception:  # pragma: no cover - fallback handled below
                schema = None

        if schema is not None:
            self._name = schema.name
            for mesh_schema in schema.meshes:
                self._load_mesh_from_schema(mesh_schema, terminal)
            self._isSet = True
        else:
            self._load_from_legacy_inputs(str(name), parameters, refine, terminal)
            self._isSet = True

        if self._isSet:
            self._finalize_integration_methods()

        if self._isSet and rank == 0:
            logging.info("Domain has been set")
            self.display()

    def _load_from_legacy_inputs(
        self, name: str, parameters: Dict[str, Any], refine: int, terminal: int
    ) -> None:
        built_in_methods = dir(scrimp.utils.mesh)

        if name in built_in_methods:
            mesh_function = getattr(scrimp.utils.mesh, name)
            gf_mesh = mesh_function(parameters, refine=refine, terminal=0)
            self._mesh.append(gf_mesh[0])
            self._dim.append(gf_mesh[1])
            self._subdomains.append(gf_mesh[2])
            self._boundaries.append(gf_mesh[3])
            self._mesh_labels.append(name)
            self._mesh_schemas.append(None)
            return

        if os.path.isfile(name):
            pathname, fileextension = os.path.splitext(name)
            basename = os.path.basename(name)[: -len(fileextension)]
            try:
                assert fileextension in [".msh", ".geo", ".py"]
                self._name = basename
                if fileextension == ".msh":
                    pass
                elif fileextension == ".geo":
                    import gmsh

                    gmsh.initialize()
                    gmsh.option.setNumber("General.Terminal", terminal)
                    gmsh.option.setNumber("General.NumThreads", 0)  # Use system default
                    gmsh.model.add(basename)
                    for key, value in parameters.items():
                        gmsh.parser.setNumber(key, value=[value])
                    gmsh.merge(name)
                    gmsh.model.geo.synchronize()
                    self._dim.append(gmsh.model.getDimension())
                    subdomains = dict()
                    for dimTags in gmsh.model.getPhysicalGroups(self._dim[-1]):
                        subdomains[
                            gmsh.model.getPhysicalName(self._dim[-1], dimTags[1])
                        ] = gmsh.model.getEntitiesForPhysicalGroup(
                            self._dim[-1], dimTags[1]
                        )[
                            0
                        ]
                    self._subdomains.append(subdomains)
                    boundaries = dict()
                    for dimTags in gmsh.model.getPhysicalGroups(self._dim[-1] - 1):
                        boundaries[
                            gmsh.model.getPhysicalName(self._dim[-1] - 1, dimTags[1])
                        ] = gmsh.model.getEntitiesForPhysicalGroup(
                            self._dim[-1] - 1, dimTags[1]
                        )[
                            0
                        ]
                    self._boundaries.append(boundaries)
                    gmsh.model.mesh.generate(gmsh.model.getDimension())
                    for i in range(refine):
                        gmsh.model.mesh.refine()
                    if rank == 0:
                        gmsh.write(os.path.join(outputs_path, "mesh", basename + ".msh"))
                    comm.barrier()
                    gmsh.finalize()
                    mesh = gf.Mesh(
                        "import",
                        "gmsh_with_lower_dim_elt",
                        os.path.join(outputs_path, "mesh", basename + ".msh"),
                    )
                    self._mesh.append(mesh)
                    self._mesh_labels.append(basename)
                    self._mesh_schemas.append(None)

                elif fileextension == ".py":
                    pass
                else:
                    logging.error(f"Construction of Domain fails from {name}")
                    raise NotImplementedError
            except AssertionError as err:
                logging.error(f"Construction of Domain fails from {name}")
                raise err
        else:
            logging.error(
                f"Construction of Domain fails from {name} (this file or function does not exist)"
            )
            raise NotImplementedError

    def _load_mesh_from_schema(self, mesh_schema: MeshSchema, terminal: int) -> None:
        source = mesh_schema.source
        generator = mesh_schema.generator

        if source == "builtin":
            built_in_methods = dir(scrimp.utils.mesh)
            if generator not in built_in_methods:
                raise ValueError(f"Unknown built-in geometry '{generator}'")
            mesh_function = getattr(scrimp.utils.mesh, generator)
            gf_mesh = mesh_function(
                mesh_schema.parameters,
                refine=mesh_schema.refinement.levels,
                terminal=terminal,
            )
            mesh = gf_mesh[0]
            dim = gf_mesh[1]
            subdomains = gf_mesh[2]
            boundaries = gf_mesh[3]
        elif source in {"gmsh", "file"}:
            path = Path(generator)
            if not path.exists():
                raise FileNotFoundError(f"Mesh generator file '{generator}' not found")
            pathname, fileextension = os.path.splitext(generator)
            basename = os.path.basename(generator)[: -len(fileextension)]
            if fileextension == ".msh":
                mesh = gf.Mesh("import", "gmsh", generator)
                dim = mesh.dim()
                subdomains = {}
                boundaries = {}
            elif fileextension == ".geo":
                import gmsh

                gmsh.initialize()
                gmsh.option.setNumber("General.Terminal", terminal)
                gmsh.option.setNumber("General.NumThreads", 0)
                gmsh.model.add(basename)
                for key, value in mesh_schema.parameters.items():
                    gmsh.parser.setNumber(key, value=[value])
                gmsh.merge(generator)
                gmsh.model.geo.synchronize()
                dim = gmsh.model.getDimension()
                subdomains = dict()
                for dimTags in gmsh.model.getPhysicalGroups(dim):
                    subdomains[
                        gmsh.model.getPhysicalName(dim, dimTags[1])
                    ] = gmsh.model.getEntitiesForPhysicalGroup(dim, dimTags[1])[0]
                boundaries = dict()
                for dimTags in gmsh.model.getPhysicalGroups(dim - 1):
                    boundaries[
                        gmsh.model.getPhysicalName(dim - 1, dimTags[1])
                    ] = gmsh.model.getEntitiesForPhysicalGroup(dim - 1, dimTags[1])[0]
                gmsh.model.mesh.generate(dim)
                for _ in range(mesh_schema.refinement.levels):
                    gmsh.model.mesh.refine()
                if rank == 0:
                    gmsh.write(os.path.join(outputs_path, "mesh", basename + ".msh"))
                comm.barrier()
                gmsh.finalize()
                mesh = gf.Mesh(
                    "import",
                    "gmsh_with_lower_dim_elt",
                    os.path.join(outputs_path, "mesh", basename + ".msh"),
                )
            else:
                raise NotImplementedError(
                    f"Unsupported mesh generator extension '{fileextension}'"
                )
        else:
            raise NotImplementedError(f"Unsupported mesh source '{source}'")

        # Apply optional aliasing from schema
        for region in mesh_schema.regions:
            subdomains[region.name] = region.tag
        for boundary in mesh_schema.boundaries:
            boundaries[boundary.name] = boundary.tag

        self._mesh.append(mesh)
        self._dim.append(dim)
        self._subdomains.append(subdomains)
        self._boundaries.append(boundaries)
        self._mesh_labels.append(mesh_schema.id)
        self._mesh_schemas.append(mesh_schema)

    def _integration_rule_for_mesh(
        self, dim: int, mesh_schema: Optional[MeshSchema]
    ) -> IntegrationRuleSchema:
        if mesh_schema and mesh_schema.integration:
            return mesh_schema.integration
        return IntegrationRuleSchema()

    def _finalize_integration_methods(self) -> None:
        self._int_method = []
        for mesh, dim, mesh_schema in zip(
            self._mesh, self._dim, self._mesh_schemas
        ):
            integration_schema = self._integration_rule_for_mesh(dim, mesh_schema)
            try:
                rule = integration_schema.as_getfem(dim)
                self._int_method.append(gf.MeshIm(mesh, gf.Integ(rule)))
            except Exception as exc:
                self._int_method.append(None)
                self._isSet = False
                if rank == 0:
                    logging.warning(
                        f"Integration method has to be set manually on mesh of dimension {dim}: {exc}"
                    )
        # Ensure boolean consistent if all integration created successfully
        if None in self._int_method:
            self._int_method = [m for m in self._int_method if m is not None]


    def set_mim_auto(self):
        """Define the integration method to a default choice"""

        self._finalize_integration_methods()

    def display(self):
        """A method giving infos about the domain"""

        try:
            assert self._isSet
        except AssertionError as err:
            logging.error("Domain has not been set yet")
            raise err

        if rank == 0:
            logging.info(f"Domain is set and contains {len(self._mesh)} mesh(es):")
            for k in range(len(self._mesh)):
                logging.info(f"=== on mesh {k} of dim {self._dim[k]}")
                logging.info(f"* Subdomains are: {self._subdomains[k]}")
                logging.info(f"* Boundaries are: {self._boundaries[k]}")

                self._mesh[k].display()  # GetFEM infos

    def get_name(self) -> str:
        """This function get the name of the domain.

        Returns:
            str: name of the domain.
        """

        return self._name

    def get_mesh(self) -> list:
        """This function gets the list of mesh for the domain.

        Returns:
            list: list of mesh for the domain
        """

        return self._mesh

    def get_mesh_labels(self) -> List[str]:
        """Return the labels used in the schema for each mesh."""

        return self._mesh_labels

    def get_schema(self) -> Optional[DomainSchema]:
        """Return the underlying domain schema if available."""

        if not self._mesh_schemas:
            return None
        if any(schema is None for schema in self._mesh_schemas):
            return None
        return DomainSchema(name=self._name, meshes=[schema for schema in self._mesh_schemas if schema is not None])

    def get_dim(self) -> list:
        """This function gets the list of dimensions for the domain.

        Returns:
            list: list of dimensions for the domain
        """

        return self._dim

    def get_subdomains(self) -> list:
        """This function gets the list of subdomains for the domain.

        Returns:
            list: list of the subdomains for the domain
        """

        return self._subdomains

    def get_boundaries(self) -> list:
        """This function gets the list of the bounderies for the domain.

        Returns:
            list: list of the boundaries for the domain
        """

        return self._boundaries

    def get_isSet(self) -> bool:
        """This function gets the boolean vale indicating wether if a Mesh has been set for the domain or not.

        Returns:
            bool: boolean indicating if a Mesh has been set
        """

        return self._isSet
