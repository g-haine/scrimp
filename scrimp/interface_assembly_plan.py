# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Read-only, fail-closed placement plan for interface weak forms.

This is an *audit*, not a GetFEM expression compiler or assembly policy.
The owner of the residual equation (Test_* variable), not the source
variable, determines the side of a cross-subdomain feedback brick.
Unknown GWFL macros/aliases and non-test expressions cannot be certified
by this restricted identifier scanner; they must be handled by a future
GWFL-aware compiler rather than blindly deployed to production.
"""

import re

from scrimp.interface_ownership import InterfaceOwnershipError

_TOKEN = re.compile(r"[A-Za-z_][A-Za-z_0-9]*")


def plan_brick_sides(interface, subdomains, port_variables, bulk_owners,
                     port_owners, bricks):
    """Plan logical-interface brick placement without mutating any object.

    port_variables: mapping of interface port -> (flow, effort)
    bulk_owners: mapping of bulk field -> subdomain region
    port_owners: inferred mapping of interface port -> subdomain region
    bricks: iterable of (brick_name, regions, form)

    Every interface brick must have exactly one recognized test variable,
    or multiple test variables with exactly one common owning subdomain.
    A no-test interface form is deliberately rejected in this stage,
    including prescribed sources and unknown macros.
    """
    interface = int(interface)
    allowed = {int(r) for r in subdomains}
    if len(allowed) != 2:
        raise InterfaceOwnershipError("Expected two distinct subdomains")
    owners = {}
    for var, region in bulk_owners.items():
        owners[var] = int(region)
    for name, variables in port_variables.items():
        if name not in port_owners:
            raise InterfaceOwnershipError(
                f"Interface port {name!r} has no certified owner"
            )
        side = int(port_owners[name])
        if side not in allowed:
            raise InterfaceOwnershipError(
                f"Port {name!r} has non-adjacent owner {side}"
            )
        for variable in variables:
            if not variable:
                continue
            if variable in owners and owners[variable] != side:
                raise InterfaceOwnershipError(
                    f"Variable {variable!r} has conflicting owners"
                )
            owners[variable] = side

    result = []
    seen = set()
    for brick_name, regions, form in bricks:
        if interface not in tuple(int(region) for region in regions):
            continue
        if brick_name in seen:
            raise InterfaceOwnershipError(f"Duplicate brick {brick_name!r}")
        seen.add(brick_name)
        tokens = set(_TOKEN.findall(form))
        tested = sorted(
            variable for variable in owners
            if "Test_" + variable in tokens
        )
        if not tested:
            raise InterfaceOwnershipError(
                f"Interface brick {brick_name!r} has no resolvable "
                f"test variable; cannot infer its assembly side"
            )
        sides = {owners[v] for v in tested}
        if len(sides) != 1:
            raise InterfaceOwnershipError(
                f"Interface brick {brick_name!r} tests variables from "
                f"multiple subdomains: {tested}"
            )
        side = next(iter(sides))
        if side not in allowed:
            raise InterfaceOwnershipError(
                f"Interface brick {brick_name!r} targets non-adjacent {side}"
            )
        result.append({
            "brick": brick_name,
            "logical_region": interface,
            "owner": side,
            "tested_variables": tested,
            "form": form,
        })
    return result
