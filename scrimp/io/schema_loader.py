"""Utilities to load SCRIMP configuration schemas from YAML or JSON files."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Literal, Optional, Tuple, Union

import yaml
from pydantic import BaseModel, Field, field_validator, model_validator


class IntegrationRuleSchema(BaseModel):
    """Describe an integration rule for a mesh."""

    name: str = Field(default="default", description="Identifier of the rule")
    expression: Optional[str] = Field(
        default=None,
        description="Explicit GetFEM integration expression (takes precedence).",
    )
    family: Literal["auto", "gauss", "triangle", "tetrahedron", "segment"] = Field(
        default="auto",
        description="Family of integration rules used when no explicit expression is provided.",
    )
    order: Optional[int] = Field(
        default=None,
        description="Order of the integration rule (if applicable).",
    )

    def as_getfem(self, dim: int) -> str:
        """Return a GetFEM integration string."""

        if self.expression:
            return self.expression

        if self.family == "auto":
            return self._auto_by_dim(dim)

        family = self.family
        order = self.order
        if family == "gauss":
            order = order or 5
            return f"IM_GAUSS1D({order})"
        if family == "triangle":
            order = order or 7
            return f"IM_TRIANGLE({order})"
        if family == "tetrahedron":
            order = order or 8
            return f"IM_TETRAHEDRON({order})"
        if family == "segment":
            order = order or 5
            return f"IM_GAUSS1D({order})"

        raise ValueError(f"Unsupported integration family: {self.family}")

    @staticmethod
    def _auto_by_dim(dim: int) -> str:
        if dim == 1:
            return "IM_GAUSS1D(5)"
        if dim == 2:
            return "IM_TRIANGLE(7)"
        if dim == 3:
            return "IM_TETRAHEDRON(8)"
        raise ValueError(f"Unsupported mesh dimension: {dim}")


class RefinementSchema(BaseModel):
    """Describe refinement strategy for a mesh."""

    levels: int = Field(default=0, ge=0)
    strategy: Literal["uniform", "auto"] = Field(default="auto")


class RegionSchema(BaseModel):
    """Describe a region (subdomain or boundary)."""

    name: str
    tag: int
    dimension: Optional[int] = None


class MeshSchema(BaseModel):
    """Describe a mesh source and annotations."""

    id: str = Field(description="Unique identifier of the mesh")
    source: Literal["builtin", "gmsh", "file"] = Field(default="builtin")
    generator: str = Field(
        description="Built-in geometry name or path to a gmsh/file generator",
    )
    parameters: Dict[str, Any] = Field(default_factory=dict)
    refinement: RefinementSchema = Field(default_factory=RefinementSchema)
    integration: Optional[IntegrationRuleSchema] = None
    regions: List[RegionSchema] = Field(default_factory=list)
    boundaries: List[RegionSchema] = Field(default_factory=list)
    description: Optional[str] = None

    @field_validator("generator")
    def strip_generator(cls, value: str) -> str:
        return value.strip()


class DomainSchema(BaseModel):
    """Top level description of the geometric domain."""

    name: str
    meshes: List[MeshSchema]

    @field_validator("meshes")
    def ensure_unique_ids(cls, meshes: List[MeshSchema]) -> List[MeshSchema]:
        ids = [m.id for m in meshes]
        if len(ids) != len(set(ids)):
            raise ValueError("Mesh ids must be unique within a domain schema")
        return meshes


ValueType = Literal["scalar", "vector", "tensor"]


class FEMFieldSchema(BaseModel):
    """Describe the discretisation of a field."""

    name: str
    order: int = 1
    family: Literal["CG", "DG", "RT", "BDM", "custom"] = "CG"
    expression: Optional[str] = Field(
        default=None,
        description="Explicit GetFEM FEM expression when family is 'custom'",
    )
    value_type: ValueType = "scalar"
    components: Optional[int] = Field(
        default=None,
        description="Number of components for vector/tensor fields",
    )
    shape: Optional[Tuple[int, int]] = Field(
        default=None,
        description="Shape for tensor-valued fields (rows, columns)",
    )
    mesh: Optional[str] = Field(
        default=None,
        description="Mesh identifier on which the FEM is defined",
    )
    description: Optional[str] = None

    @model_validator(mode="after")
    def check_expression(cls, values: "FEMFieldSchema") -> "FEMFieldSchema":
        family = values.family
        expression = values.expression
        if family == "custom" and not expression:
            raise ValueError("Custom FEM requires an explicit expression")
        return values

    @field_validator("shape")
    def ensure_positive_shape(
        cls, shape: Optional[Tuple[int, int]]
    ) -> Optional[Tuple[int, int]]:
        if shape is None:
            return shape
        rows, cols = shape
        if rows <= 0 or cols <= 0:
            raise ValueError("Tensor shape dimensions must be positive")
        return shape

    def value_size(self, mesh_dim: Optional[int] = None) -> int:
        """Infer the number of components for the discretisation."""

        if self.value_type == "scalar":
            return 1
        if self.value_type == "vector":
            if self.components is not None:
                return self.components
            if mesh_dim is None:
                raise ValueError(
                    "Vector FEM requires 'components' or mesh dimension context"
                )
            return mesh_dim
        if self.value_type == "tensor":
            if self.shape is not None:
                return self.shape[0] * self.shape[1]
            if self.components is not None:
                return self.components
            if mesh_dim is None:
                raise ValueError(
                    "Tensor FEM requires 'components' or mesh dimension context"
                )
            return mesh_dim * mesh_dim
        raise ValueError(f"Unsupported value_type {self.value_type}")


class FEMSchema(BaseModel):
    """Container of field discretisations."""

    fields: List[FEMFieldSchema] = Field(default_factory=list)


class MaterialParameterSchema(BaseModel):
    """Describe a material parameterisation."""

    name: str
    value: Union[float, int, List[float], Dict[str, Any]]
    units: Optional[str] = None
    description: Optional[str] = None


class SimulationSchema(BaseModel):
    """Root schema bundling domain, FEM and material definitions."""

    domain: DomainSchema
    fem: FEMSchema = Field(default_factory=FEMSchema)
    materials: Dict[str, MaterialParameterSchema] = Field(default_factory=dict)


def _load_yaml_or_json(path: Path) -> Dict[str, Any]:
    text = path.read_text()
    if path.suffix.lower() in {".yaml", ".yml"}:
        return yaml.safe_load(text)
    if path.suffix.lower() == ".json":
        return json.loads(text)
    raise ValueError(f"Unsupported schema format: {path.suffix}")


def load_schema(source: Union[str, Path, Dict[str, Any]]) -> SimulationSchema:
    """Load a :class:`SimulationSchema` from a file path or dictionary."""

    if isinstance(source, (str, Path)):
        data = _load_yaml_or_json(Path(source))
    else:
        data = source
    return SimulationSchema.model_validate(data)


def load_domain_schema(source: Union[str, Path, Dict[str, Any], DomainSchema]) -> DomainSchema:
    if isinstance(source, DomainSchema):
        return source
    if isinstance(source, (str, Path)):
        data = _load_yaml_or_json(Path(source))
    else:
        data = source
    if "domain" in data:
        return DomainSchema.model_validate(data["domain"])
    return DomainSchema.model_validate(data)


def load_fem_schema(source: Union[str, Path, Dict[str, Any], FEMSchema]) -> FEMSchema:
    if isinstance(source, FEMSchema):
        return source
    if isinstance(source, (str, Path)):
        data = _load_yaml_or_json(Path(source))
    else:
        data = source
    if "fem" in data:
        return FEMSchema.model_validate(data["fem"])
    return FEMSchema.model_validate(data)


# ---------------------------------------------------------------------------
# Built-in templates
# ---------------------------------------------------------------------------


def _builtin_domain_templates() -> Dict[str, DomainSchema]:
    return {
        "interval": DomainSchema(
            name="Interval",
            meshes=[
                MeshSchema(
                    id="omega",
                    source="builtin",
                    generator="Interval",
                    parameters={"L": 1.0, "h": 0.05},
                    description="One-dimensional interval with uniform spacing",
                )
            ],
        ),
        "rectangle": DomainSchema(
            name="Rectangle",
            meshes=[
                MeshSchema(
                    id="omega",
                    source="builtin",
                    generator="Rectangle",
                    parameters={"L": 2.0, "l": 1.0, "h": 0.1},
                    description="Two-dimensional rectangle domain",
                )
            ],
        ),
        "disk": DomainSchema(
            name="Disk",
            meshes=[
                MeshSchema(
                    id="omega",
                    source="builtin",
                    generator="Disk",
                    parameters={"R": 1.0, "h": 0.1},
                    description="Circular domain for axisymmetric problems",
                )
            ],
        ),
    }


def _builtin_material_templates() -> Dict[str, MaterialParameterSchema]:
    return {
        "heat_homogeneous": MaterialParameterSchema(
            name="thermal_conductivity",
            value=205.0,
            units="W/(m*K)",
            description="Default aluminium conductivity",
        ),
        "wave_density": MaterialParameterSchema(
            name="density",
            value=7800.0,
            units="kg/m^3",
            description="Steel density for wave propagation",
        ),
    }


def _builtin_fem_templates() -> Dict[str, FEMSchema]:
    return {
        "scalar_lagrange": FEMSchema(
            fields=[
                FEMFieldSchema(
                    name="u",
                    family="CG",
                    order=1,
                    value_type="scalar",
                )
            ]
        ),
        "vector_lagrange": FEMSchema(
            fields=[
                FEMFieldSchema(
                    name="v",
                    family="CG",
                    order=1,
                    value_type="vector",
                )
            ]
        ),
    }


DEFAULT_DOMAIN_LIBRARY = _builtin_domain_templates()
DEFAULT_MATERIAL_LIBRARY = _builtin_material_templates()
DEFAULT_FEM_LIBRARY = _builtin_fem_templates()


def list_templates(kind: str) -> Iterable[str]:
    if kind == "domain":
        return DEFAULT_DOMAIN_LIBRARY.keys()
    if kind == "material":
        return DEFAULT_MATERIAL_LIBRARY.keys()
    if kind == "fem":
        return DEFAULT_FEM_LIBRARY.keys()
    raise ValueError(f"Unknown template kind: {kind}")


# ---------------------------------------------------------------------------
# CLI helpers
# ---------------------------------------------------------------------------


def _cmd_list(args: argparse.Namespace) -> None:
    items = list_templates(args.kind)
    for item in items:
        print(item)


def _cmd_generate(args: argparse.Namespace) -> None:
    domain = DEFAULT_DOMAIN_LIBRARY.get(args.domain)
    if domain is None:
        raise SystemExit(f"Unknown domain template '{args.domain}'")

    fem_schema = DEFAULT_FEM_LIBRARY.get(args.fem, FEMSchema())
    materials: Dict[str, MaterialParameterSchema] = {}
    if args.material is not None:
        material = DEFAULT_MATERIAL_LIBRARY.get(args.material)
        if material is None:
            raise SystemExit(f"Unknown material template '{args.material}'")
        materials[material.name] = material

    schema = SimulationSchema(domain=domain, fem=fem_schema, materials=materials)
    payload = schema.dict()
    text = yaml.safe_dump(payload, sort_keys=False)
    if args.output:
        Path(args.output).write_text(text)
    else:
        print(text)


def _cmd_validate(args: argparse.Namespace) -> None:
    schema = load_schema(args.path)
    print(schema.json(indent=2))


def main(argv: Optional[List[str]] = None) -> None:
    parser = argparse.ArgumentParser(description="SCRIMP schema helper")
    sub = parser.add_subparsers(dest="command", required=True)

    list_cmd = sub.add_parser("list", help="List available templates")
    list_cmd.add_argument("kind", choices=["domain", "material", "fem"])
    list_cmd.set_defaults(func=_cmd_list)

    gen_cmd = sub.add_parser("generate", help="Generate a schema from templates")
    gen_cmd.add_argument("domain", help="Domain template name")
    gen_cmd.add_argument("--fem", default="scalar_lagrange", help="FEM template name")
    gen_cmd.add_argument(
        "--material",
        default=None,
        help="Material template name (optional)",
    )
    gen_cmd.add_argument(
        "--output",
        default=None,
        help="Optional output path (defaults to stdout)",
    )
    gen_cmd.set_defaults(func=_cmd_generate)

    val_cmd = sub.add_parser("validate", help="Validate an existing schema")
    val_cmd.add_argument("path", help="Path to schema file")
    val_cmd.set_defaults(func=_cmd_validate)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
