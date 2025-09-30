Builder quickstart
==================

The new builder API lets you register domain, states and ports in a fluent fashion::

    from scrimp import CoState, DPHS, Domain, State

    domain = Domain("Interval", {"L": 1.0, "h": 0.2})
    dphs = DPHS()

    state = State("x", "displacement", "scalar-field")
    costate = CoState("p", "momentum", state)

    dphs.builder.with_domain(domain).add_state(state).add_costate(costate)

    print(list(dphs.ports))

See :mod:`examples.builder_quickstart` for a complete script.
