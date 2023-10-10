#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import scrimp as S

## Parameters
# Space
L = 1
h = 0.002
FE_Phi_order = 1
FE_Psi_order = 2
FE_Flux_order = 1

# Time scheme
t_0 = 0.0
t_f = 2.0
dt = h/10.
dt_save = 0.01
ts_adapt_dt_min = h/1000.
ts_type = "bdf"
bdf_order = 4
ksp_type = "preonly"
pc_type = "lu"
init_step = True

# Physics
Eta = 2
epsilon = 0.002
w = 0.000005
Delta_s = 1

# Controls
U_L = "0."
U_R = "0."

# Initial values
C = 36
Phi_init = f"1/(1+np.exp(-{C}*(x-1/2)))-1/np.exp(-{C}*(0-1/2))"
Psi_init = f"{C}*np.exp(-{C}*(x-1/2))/((np.exp(-{C}*(x-1/2))+1)*(np.exp(-{C}*(x-1/2))+1))"

# Real distributed pHs initialisation
AC = S.DPHS("real")

# Define the domain...
AC.set_domain(S.Domain('Interval', {'L': L, 'h': h}))

# Define the states...
States = [
        S.State('Phi', 'Phase field', 'scalar-field'),
        S.State('Psi', 'grad(Phi)', 'scalar-field'),
    ]

# Loop on the list States
for state in States:
    # Add the state
    AC.add_state(state)
    # And its co-state
    AC.add_costate(S.CoState('e_'+state.get_name(), 'co-state '+state.get_name(), state))

# Define the algebraic port
AC.add_port(S.Port('Phase field flux', 'F_AC', 'E_AC', 'scalar-field', dissipative=False))

# Define the control port
AC.add_control_port(S.Control_Port('Boundary control (left)', 'U_L', 'E_AC', 'Y_L', 'e_Psi', 'scalar-field', region=10, position='effort'))
AC.add_control_port(S.Control_Port('Boundary control (right)', 'U_R', 'E_AC', 'Y_R', 'e_Psi', 'scalar-field', region=11, position='effort'))

# Define the finite element spaces...
FE_Phi = S.FEM('Phi', FE_Phi_order)
FE_Psi = S.FEM('Psi', FE_Psi_order)
FE_Flux = S.FEM('Phase field flux', FE_Flux_order)
FE_left = S.FEM('Boundary control (left)', 1)
FE_right = S.FEM('Boundary control (right)', 1)
# ... and plug them in the pHs
AC.add_FEM(FE_Phi)
AC.add_FEM(FE_Psi)
AC.add_FEM(FE_Flux)
AC.add_FEM(FE_left)
AC.add_FEM(FE_right)

# Add parameters
AC.add_parameter(S.Parameter('Eta', 'Interface mobility', 'scalar-field', f'{Eta}', 'Phase field flux'))
AC.add_parameter(S.Parameter('epsilon', 'related to the thickness of the interface', 'scalar-field', f'{epsilon}', 'Psi'))
AC.add_parameter(S.Parameter('w', ' tuning of the hallow between the two local maximum', 'scalar-field', f'{w}', 'Phase field flux'))
AC.add_parameter(S.Parameter('Delta_s', 'diff between liquid-solid ', 'scalar-field', f'{Delta_s}', 'Phase field flux'))

# Define polynomials
p_i = lambda phi: f"6*{phi}*{phi}*{phi}*{phi}*{phi} - 15*{phi}*{phi}*{phi}*{phi} + 10*{phi}*{phi}*{phi}"
expr_p_i = p_i("Phi")

p_w = lambda phi: f" -{phi}*{phi}*{phi}*{phi} + 2*{phi}*{phi}*{phi} - {phi}*{phi}"
expr_p_w = p_w("Phi")

AC.add_brick(S.Brick('M_Phi', 'Phi*Test_Phi', [1], dt=True, position='flow'))
AC.add_brick(S.Brick('M_Psi', 'Psi*Test_Psi', [1], dt=True, position='flow'))
AC.add_brick(S.Brick('M_F', 'F_AC.Test_F_AC', [1], position='flow'))
AC.add_brick(S.Brick('M_Y_L', 'Y_L*Test_Y_L', [10], position='flow'))
AC.add_brick(S.Brick('M_Y_R', 'Y_R*Test_Y_R', [11], position='flow'))

AC.add_brick(S.Brick('M_E', '-E_AC*Test_Phi', [1], position='effort'))

AC.add_brick(S.Brick('D', 'E_AC*Div(Test_Psi)', [1], position='effort'))

AC.add_brick(S.Brick('B_L', 'U_L*Test_Psi', [10], position='effort'))
AC.add_brick(S.Brick('B_R', '-U_R*Test_Psi', [11], position='effort'))

AC.add_brick(S.Brick('M_I', 'e_Phi*Test_F_AC', [1], position='effort'))
AC.add_brick(S.Brick('-D^T', '-Div(e_Psi)*Test_F_AC', [1], position='effort'))

AC.add_brick(S.Brick('C_L', '-e_Psi*Test_Y_L', [10], position='effort'))
AC.add_brick(S.Brick('C_R', 'e_Psi*Test_Y_R', [11], position='effort'))

AC.add_brick(S.Brick('-M_E_AC', '-E_AC.Test_E_AC', [1]))
AC.add_brick(S.Brick('CR_AC', 'F_AC/Eta.Test_E_AC', [1]))

AC.add_brick(S.Brick('-M_e_Phi', '-e_Phi.Test_e_Phi', [1]))
AC.add_brick(S.Brick('CR_e_Phi_1', f"Delta_s*Diff({expr_p_i}, Phi, Test_e_Phi)", [1], linear=False, explicit=False))
AC.add_brick(S.Brick('CR_e_Phi_2', f"w*Diff({expr_p_w}, Phi, Test_e_Phi)", [1], linear=False, explicit=False))

AC.add_brick(S.Brick('-M_e_Psi', '-e_Psi*Test_e_Psi', [1]))
AC.add_brick(S.Brick('CR_e_Psi', 'epsilon*Psi*Test_e_Psi', [1]))

# Set the control functions
AC.set_control('Boundary control (left)', U_L)
AC.set_control('Boundary control (right)', U_R)

# Set the initial data
AC.set_initial_value('Phi', Phi_init)
AC.set_initial_value('Psi', Psi_init)

# Solve
AC.set_time_scheme(t_0=t_0, t_f=t_f, dt=dt, ts_adapt_dt_min=ts_adapt_dt_min,
                   ts_type=ts_type, bdf_order=bdf_order,
                   ksp_type=ksp_type, pc_type=pc_type, init_step=init_step)
AC.solve()

# Hamiltonian
AC.hamiltonian.add_term(S.Term('Entropy', f"Delta_s*{expr_p_i} + w*{expr_p_w}", [1]))
AC.hamiltonian.add_term(S.Term('Interface energy', "-0.5*epsilon*Psi*Psi", [1]))
AC.hamiltonian.set_name('Hamiltonian \overline{S}')

AC.plot_Hamiltonian(with_powers=True)

import matplotlib.pyplot as plt  
import numpy as np  
# dofs_Phi = AC.gf_model.interval_of_variable("Phi")  

t = np.array(AC.solution['t'])  
x = np.linspace(0,L,int(1/h)+1)
Phi = AC.get_solution("Phi")

# np.array([AC.solution['z'][k].array[dofs_Phi[0]:dofs_Phi[0]+dofs_Phi[1]] for k in range(len(t))])  

# sol = np.array([AC.solution['z'][k].array[dofs_Phi[0]:dofs_Phi[0]+dofs_Phi[1]] for k in [0,int(len(t)/10),int(1*len(t)/2),int(3*len(t)/4)]])  

sol = np.array([Phi[k] for k in [0,int(len(t)/10),int(1*len(t)/2),int(3*len(t)/4)]])

plt.clf()  
plt.plot(x,np.transpose(sol))
plt.xlabel('Linear space x')  
plt.ylabel("Phi")  
plt.title("Phi over Linear space")  
plt.legend([f"t={t[tk]:8g}" for tk in [0,int(len(t)/10),int(1*len(t)/2),int(3*len(t)/4)]])
plt.show()

# Export solutions for ParaView  
# AC.export_to_pv('Phi')  
# AC.export_to_pv('E_AC')  

## Show the evolution of the phase over the domain
# prove to be constant for Cahn-Hilliard, see our GSI 2023
# must grow here for a solidification process

int_Phi = AC.get_quantity("Phi")

relative = [int_Phi[k]/int_Phi[0] for k in range(len(int_Phi))]

fig = plt.figure(figsize=[8,5])
ax1 = fig.add_subplot(211)
ax1.set_title("$\int_\Omega \phi$, evolution over time")
ax1.plot(t, int_Phi)
ax1.grid(axis='both')
ax1.set_ylabel('$\int \phi$')

ax2 = fig.add_subplot(212,sharex=ax1)
ax2.plot(t, relative)
ax2.grid(axis='both')
ax2.set_xlabel('Time t')
ax2.set_ylabel('Relatively to $\phi(0)$')
plt.show()
