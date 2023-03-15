from src.dphs import DPHS
from src.domain import Domain
from src.dphs import Term
from src.state import State
from src.costate import CoState
from src.port import Parameter, Port
from src.brick import Brick

wave = DPHS("real")

wave.set_domain(Domain("Interval", {"L": 1.0, "h": 0.01}))

state = State("q", "Strain", "scalar-field")
wave.add_state(state)
wave.add_costate(CoState("e_q", "Stress", state))
state = State("p", "Linear momentum", "scalar-field")
wave.add_state(state)
wave.add_costate(CoState("e_p", "velocity", state))

wave.add_control_port(
    "Boundary control (left)",
    "U_L",
    "Velocity",
    "Y_L",
    "Normal force",
    "scalar-field",
    region=10,
)
wave.add_control_port(
    "Boundary control (right)",
    "U_R",
    "Velocity",
    "Y_R",
    "Normal force",
    "scalar-field",
    region=11,
)

wave.add_FEM("q", 2)
wave.add_FEM("p", 1)
wave.add_FEM("Boundary control (left)", 1)
wave.add_FEM("Boundary control (right)", 1)

wave.add_brick(Brick("M_q", "q.Test_q", [1], dt=True, position="flow"))
wave.add_brick(Brick("M_p", "p*Test_p", [1], dt=True, position="flow"))
wave.add_brick(Brick("M_Y_L", "Y_L*Test_Y_L", [10], position="flow"))
wave.add_brick(Brick("M_Y_R", "Y_R*Test_Y_R", [11], position="flow"))

wave.add_brick(Brick("D", "Grad(e_p).Test_q", [1], position="effort"))

wave.add_brick(Brick("-D^T", "-e_q.Grad(Test_p)", [1], position="effort"))
wave.add_brick(Brick("B_L", "U_L*Test_p", [10], position="effort"))
wave.add_brick(Brick("B_R", "U_R*Test_p", [11], position="effort"))

wave.add_brick(Brick("-B_L^T", "-e_p*Test_Y_L", [10], position="effort"))
wave.add_brick(Brick("-B_R^T", "-e_p*Test_Y_R", [11], position="effort"))

wave.add_parameter(Parameter("T", "Young's modulus", "scalar-field", "1", "q"))
wave.add_parameter(Parameter("rho", "Mass density", "scalar-field", "1 + x*(1-x)", "p"))

wave.add_brick(Brick("-M_e_q", "-e_q.Test_e_q", [1]))
wave.add_brick(Brick("CR_q", "q.T*Test_e_q", [1]))

wave.add_brick(Brick("-M_e_p", "-e_p*Test_e_p", [1]))
wave.add_brick(Brick("CR_p", "p/rho*Test_e_p", [1]))

wave.set_control("Boundary control (left)", "sin(2*pi*t)")
wave.set_control("Boundary control (right)", "0.")

wave.set_initial_value("q", "10.")
wave.set_initial_value("p", "np.exp(-50.*(x-0.5)*(x-0.5))")

wave.solve()

wave.hamiltonian.add_term(Term("Kinetic energy", "0.5*p*p/rho", [1]))
wave.hamiltonian.add_term(Term("Potential energy", "0.5*q*T*q", [1]))

wave.plot_Hamiltonian()
