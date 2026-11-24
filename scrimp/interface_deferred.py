# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
#
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""Production, conservative GetFEM compilation of conforming interfaces.

Enabled automatically when a two-sided conforming internal interface
has an unambiguous geometric label and weak-form owner. The original
public regions are immutable; register all FEs/bricks in their
declared order after resolving each parent-convex outward normal.

This compiler FAILS CLOSED on ambiguous, nonconforming, multi-labeled
or unsupported forms; it is not a universal H(div)/H(curl) trace
compiler. Two concentric Heat-Wave region permutations were qualified
with real K/DAE/CN/energy MPI1/2 and physical VTK. Other families and
geometries are not automatically certified by that evidence.
SCRIMP_AUTO_INTERFACE=0 opts out; legacy PILOT env flag is accepted.
"""

from scrimp.interface_discovery import discover_conforming_interfaces


class InterfaceCompilationError(ValueError):
    """The declared model cannot safely be compiled with oriented faces."""


class DeferredInterfaceCompiler:
    def __init__(self, domain):
        self.interfaces = discover_conforming_interfaces(domain)
        self.fems = []
        self.bricks = []
        self.replaying = False
        self.completed = False
        self.report = None

    def record_fem(self, fem):
        if self.completed:
            raise InterfaceCompilationError(
                "Cannot add FEM after first automatic interface compilation"
            )
        self.fems.append(fem)

    def record_brick(self, brick):
        if self.completed:
            raise InterfaceCompilationError(
                "Cannot add brick after first automatic interface compilation"
            )
        self.bricks.append(brick)

    def compile(self, system):
        """Validate the full model *before* registering any deferred form."""
        if self.completed:
            return self.report
        if self.replaying:
            raise InterfaceCompilationError("Recursive interface compilation")
        if not self.interfaces:
            raise InterfaceCompilationError(
                "No internal interface topology was discovered"
            )

        planned = {}
        for (mesh_id, logical), sides in sorted(self.interfaces.items()):
            audit = system.audit_interface_brick_plan(
                logical, sides.subdomains, mesh_id=mesh_id,
            )
            owners = audit["port_owners"]
            bricks = audit["brick_plan"]
            # Brick names are globally unique in DPHS.bricks. Every
            # physical interface form must have one tested-variable owner.
            if not bricks:
                raise InterfaceCompilationError(
                    f"Interface {logical} has no resolvable weak forms"
                )
            planned[(mesh_id, logical)] = (sides, owners, bricks)

        # Commit only after a successful, fail-closed preflight. Reserve
        # new region labels above all existing physical group labels.
        regions = {}
        manifests = []
        for (mesh_id, logical), (sides, owners, bricks) in sorted(planned.items()):
            mesh = system.domain.get_mesh()[mesh_id]
            taken = set(int(x) for x in mesh.regions())
            for previous in regions:
                if previous[0] == mesh_id:
                    taken.update(regions[previous].values())
            next_id = max(taken | {0}) + 1
            assignment = {}
            for bulk in sorted(sides.subdomains):
                while next_id in taken:
                    next_id += 1
                mesh.set_region(next_id, sides.for_subdomain(bulk))
                assignment[bulk] = next_id
                taken.add(next_id)
                next_id += 1
            regions[(mesh_id, logical)] = assignment
            manifests.append({
                "mesh_id": mesh_id, "logical_region": logical,
                "subdomains": list(sides.subdomains),
                "facet_count": len(sides.facet_keys),
                "side_regions": dict(assignment),
                "port_owners": dict(owners),
                "brick_owners": {
                    item["brick"]: int(item["owner"]) for item in bricks
                },
            })

        for (mesh_id, logical), (sides, owners, bricks) in planned.items():
            assignment = regions[(mesh_id, logical)]
            for port_name, bulk in owners.items():
                port = system.ports[port_name]
                if port.get_mesh_id() != mesh_id or port.get_region() != logical:
                    raise InterfaceCompilationError(
                        f"Port {port_name} declaration changed during compilation"
                    )
                port.set_assembly_region(assignment[bulk])
            for entry in bricks:
                brick = system.bricks[entry["brick"]]
                if brick.get_mesh_id() != mesh_id:
                    raise InterfaceCompilationError(
                        f"Brick {entry['brick']} mesh changed"
                    )
                original = tuple(brick.get_regions())
                if original != (logical,):
                    raise InterfaceCompilationError(
                        f"Brick {entry['brick']} spans multiple regions; "
                        "not yet supported by automatic interface compiler"
                    )
                brick.set_assembly_regions([assignment[entry["owner"]]])

        self.replaying = True
        try:
            # Replay in the *original* FEM/Brick declaration order:
            # native variable intervals and the full sparse K matrix
            # should remain directly comparable to the reference.
            for fem in self.fems:
                system.add_FEM(fem)
            # Parameters declared before an FEM was registered are now
            # ready for initialization in their original port's basis.
            for port in system.ports.values():
                if port.get_isSet():
                    for parameter in port.get_parameters():
                        system.init_parameter(parameter.get_name(), port.get_name())
            for brick in self.bricks:
                system.add_brick(brick)
        finally:
            self.replaying = False

        self.completed = True
        self.report = {
            "schema": "scrimp-v130-deferred-auto-internal-interface-v1",
            "status": "COMPILED_CONFORMING_INTERFACE",
            "interfaces": manifests,
            "fem_count": len(self.fems),
            "brick_count": len(self.bricks),
            "declaration_order_preserved": True,
            "original_logical_region_labels_preserved": True,
            # These flags distinguish the generic compiler's
            # unqualified configurations from its proven Heat-Wave
            # reference; successful construction alone never proves a
            # model's physical power balance or MPI2 equivalence.
            "physically_qualified": False,
            "MPI2_qualified": False,
            "reference_qualification": (
                "2D concentric Heat-Wave, heat/wave=1/2 and 2/1; "
                "real strict DAE, CN, pH energy, MPI1/2 and rank0 VTK"
            ),
            "production_default": True,
        }
        return self.report
