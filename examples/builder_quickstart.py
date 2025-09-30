"""Minimal example showcasing the declarative builder workflow."""

from scrimp import CoState, DPHS, Domain, State


def build_single_state_model() -> DPHS:
    domain = Domain("Interval", {"L": 1.0, "h": 0.2})
    dphs = DPHS()

    state = State("x", "displacement", "scalar-field", mesh_id=0)
    costate = CoState("p", "momentum", state)

    dphs.builder.with_domain(domain).add_state(state).add_costate(costate)

    return dphs


if __name__ == "__main__":
    model = build_single_state_model()
    print("Registered ports:", list(model.ports))
