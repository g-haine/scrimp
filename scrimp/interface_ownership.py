# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2026 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp


"""Conservative, side-effect-free ownership audit for interface ports.

This stage does NOT change GetFEM variable registration or brick assembly.
It derives ownership only from explicit bulk-variable anchors in a weak
form. Interface-only feedback is never used to infer ownership: those
terms can intentionally couple opposite subdomains.
"""

import re


class InterfaceOwnershipError(ValueError):
    """Port-side ownership was ambiguous or contradictory."""


_TOKEN = re.compile(r"[A-Za-z_][A-Za-z_0-9]*")


def variables_in_form(form, variables):
    """Return variables mentioned as identifiers or Test_<variable>.

    This is a conservative audit, not a complete GetFEM GWFL parser.
    Variables in macros and unrecognized syntax require a separate
    formal parsing/validation stage before automatic assembly.
    """
    tokens = set(_TOKEN.findall(form))
    return {
        name for name in variables
        if name in tokens or ("Test_" + name) in tokens
    }


def infer_interface_port_owners(interface, port_variables, bulk_owners, bricks):
    """Infer side owners for interface ports from anchored weak forms.

    Parameters
    ----------
    interface : int
        Logical region ID.
    port_variables : dict[str, tuple[str, str]]
        Interface-port name -> (flow variable, effort variable).
    bulk_owners : dict[str, int]
        Bulk-variable name -> its physical volume-region ID.
    bricks : iterable[tuple[iterable[int], str]]
        Immutable (regions, original_form) declarations.

    Returns
    -------
    dict[str, int]
        Interface-port name -> unique owning bulk region.

    This function uses only bricks containing explicit bulk variables
    as anchors. Interface-only mass and cross-feedback forms do not
    propagate orientation. Every port must have exactly one anchored
    side; otherwise this method fails closed.
    """
    interface = int(interface)
    variables = set(bulk_owners)
    per_port = {}
    for name, fields in port_variables.items():
        fields = tuple(fields)
        if not fields or any(not isinstance(v, str) or not v for v in fields):
            raise InterfaceOwnershipError(f"Invalid variables for port {name!r}")
        if any(v in bulk_owners for v in fields):
            raise InterfaceOwnershipError(
                f"Interface port {name!r} reuses a bulk variable"
            )
        per_port[name] = fields
        variables.update(fields)

    owners = {name: set() for name in per_port}
    for regions, form in bricks:
        if interface not in tuple(int(r) for r in regions):
            continue
        used = variables_in_form(form, variables)
        bulk_sides = {int(bulk_owners[v]) for v in used & bulk_owners.keys()}
        if not bulk_sides:
            # Cannot decide orientation from boundary-only feedback,
            # mass terms, or prescribed terms. Do not guess.
            continue
        if len(bulk_sides) != 1:
            raise InterfaceOwnershipError(
                f"Interface {interface}: form contains bulk variables "
                f"from multiple subdomains: {form!r}"
            )
        side = next(iter(bulk_sides))
        for name, fields in per_port.items():
            if used.intersection(fields):
                owners[name].add(side)

    result = {}
    for name, sides in owners.items():
        if len(sides) != 1:
            raise InterfaceOwnershipError(
                f"Interface {interface}: port {name!r} has "
                f"{'no anchored owner' if not sides else 'conflicting owners'}"
            )
        result[name] = next(iter(sides))
    return result
