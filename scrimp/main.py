from dphs.state import State
from dphs.costate import CoState
from dphs.port import Port, Parameter
from dphs.domain import Domain


state = State("name_state", "description_state", "kind_state")
costate = CoState("name_costate", "description_costate", state)
port = Port("name_port", "flow_port", "effort_port", "kind_port", 0, True, True, 1)
parameter = Parameter(
    "name_param", "description_param", "kind_param", "expression_param", port.get_name()
)
domain = Domain("Rectangle", {"L": 2.0, "l": 1.0, "h": 0.15})
wave.set_domain("Rectangle", {"L": 2.0, "l": 1.0, "h": 0.15})

## Define the variables and their discretizations

# Add a state
wave.add_state("q", "Strain", "vector-field")
# Add its co-state
wave.add_costate("e_q", "Stress", "q")
# Add a Finite Element Method to the `port`
wave.add_FEM("q", 1, FEM="DG")
# Add a (possibly space-varying) parameter to the `port`
wave.add_parameter("T", "Young's modulus", "tensor-field", "[[5+x,x*y],[x*y,2+y]]", "q")


print(state)
print(costate)
print(costate.get_state())
print(state.get_costate())
state.set_costate(costate)
print(state.get_costate())
state.set_port(port)
