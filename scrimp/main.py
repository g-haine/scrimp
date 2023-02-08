from scrimp.dphs.state import State
from scrimp.dphs.costate import CoState


s = State("1","1","1")
c = CoState("2","2",s)
print(s)
print(c)
print(c.get_state())
print(s.get_costate())
s.set_costate(c)
print(s.get_costate())
