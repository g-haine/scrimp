# SCRIMP - Simulation and ContRol of Interactions in Multi-Physics
#
# Copyright (C) 2015-2023 ISAE-SUPAERO -- GNU GPLv3
# 
# See the LICENSE file for license information.
#
# github: https://github.com/g-haine/scrimp

"""
- file:             sandbox/INSE_vorticity_stream.py
- authors:          Ghislain Haine
- date:             27 jun. 2023
- brief:            2D Incompressible Navier-Stokes Equations 
                    in vorticity-stream function formulation
                    Lid-driven cavity at Reynolds 100
"""

def INSE():
    L = 1.0
    l = 1.0
    h = 0.1
    hmin = 0.05
    layer = 0.1
    refine = 0
    
    order_psi = 3
    order_d = 2
    order_bnd = 1
    
    t_f = 10.
    dt = 0.01
    dt_save = 0.1
    ts_adapt_dt_min = 1e-5
    ts_type = "bdf"
    bdf_order = 4
    ksp_type = "preonly"
    solver_type = "mumps"
    init_step = True
    
    rho = 1.
    Reynolds = 400.
    mu = 1./Reynolds
    
    boundaries = ["bottom", "right", "top", "left"]
    control_types = ["dt", "normal", "curl", "tangent"]
    control_expressions = [["0.", "0.5", "0.", "0."],
                           ["0.", "0.5", "0.", "0."],
                           ["0.", "0.5", "0.", "1."],
                           ["0.", "0.5", "0.", "0."]]
    
    psi_0 = "0.5"
    
    import scrimp as S
    
    inse = S.DPHS("real")
    
    inse.gf_model.add_macro("Rot", "[[0,1],[-1,0]]")
    inse.gf_model.add_macro("GradPerp(v)", "Rot*Grad(v)")
    
    inse.set_domain(S.Domain("Rectangle.geo", {"L": L, "l": l, "h": h, "hmin": hmin, "layer": layer}, refine=refine))
    
    psi = S.State("psi", "Stream function", "scalar-field")
    inse.add_state(psi)
    inse.add_costate(S.CoState("e_psi", "Stream function", psi, substituted=True))
    inse.add_FEM(S.FEM("psi", order_psi))
    
    d = S.Port("Dissipation", "e_d", "e_d", "scalar-field", substituted=True, dissipative=True)
    inse.add_port(d)
    inse.add_FEM(S.FEM("Dissipation", order_d))
    
    # Loop for the 16 control ports
    control_ports = list()
    for boundary_k in range(len(boundaries)):
        boundary = boundaries[boundary_k]
        for control_type in control_types:
            if control_type == "curl":
                position = "flow"
            else:
                position = "effort"
            control_ports.append(S.Control_Port("Control "+control_type+" "+boundary, 
                                                "U_"+control_type[0]+"_"+boundary[0], 
                                                control_type+" control on "+boundary, 
                                                "Y_"+control_type[0]+"_"+boundary[0], 
                                                control_type+" observation on "+boundary, 
                                                "scalar-field", region=int(boundary_k+10), 
                                                position=position))
    
    for control_port in control_ports:
        inse.add_control_port(control_port)
        inse.add_FEM(S.FEM(control_port.get_name(), order_bnd))
    
    inse.add_parameter(S.Parameter("rho", "Mass density", "scalar-field", f"{rho}", "psi"))
    inse.add_parameter(S.Parameter("mu", "Viscosity", "scalar-field", f"{mu}", "Dissipation"))
    
    bricks = [
            S.Brick("M_omega", "rho * Grad(psi) . Grad(Test_psi)", [1], dt=True, position="flow"),
            S.Brick("M_c", "e_d / mu * Test_e_d", [1], position="flow")
        ]
    
    # Loop for the 16 mass matrices, 16 control matrices and 16 observation matrices
    for control_port in control_ports:
        if control_port.get_position() == "effort":
            name_control = control_port.get_name_control()
            name_observation = control_port.get_name_obervation()
            bricks.append(S.Brick(
                    "M"+name_observation[1:], name_observation+"*Test_"+name_observation, 
                    [control_port.get_region()], position="flow"
                ))
            if "tangent" in control_port.get_name():
                field = "-e_d"
            else:
                field = " psi"
            bricks.append(S.Brick(
                    "B"+name_control[1:], name_control+"*Test_"+field[1:], 
                    [control_port.get_region()], position="effort"
                ))
            bricks.append(S.Brick(
                    "-B"+name_control[1:]+"^T", "-"+field+"*Test_"+name_observation, 
                    [control_port.get_region()], position="effort"
                ))
        else:
            name_control = control_port.get_name_control()
            name_observation = control_port.get_name_obervation()
            bricks.append(S.Brick(
                    "M"+name_observation[1:], name_control+"*Test_"+name_observation, 
                    [control_port.get_region()], position="flow"
                ))
            bricks.append(S.Brick(
                    "B"+name_control[1:], "-"+name_observation+"*Test_psi", 
                    [control_port.get_region()], position="effort"
                ))
            bricks.append(S.Brick(
                    "-B"+name_control[1:]+"^T", "psi*Test_"+name_observation, 
                    [control_port.get_region()], position="effort"
                ))
    
    bricks.append(S.Brick(
            "J(omega)", "e_d/mu * GradPerp(psi) . Grad(Test_psi)", [1], 
            linear=False, explicit=False, position="effort"
        ))
    
    bricks.append(S.Brick(
            "-D", "- Grad(e_d) . Grad(Test_psi)", [1], position="effort"
        ))
    
    bricks.append(S.Brick(
            "D^T", "Grad(psi) . Grad(Test_e_d)", [1], position="effort"
        ))
    
    for brick in bricks:
        inse.add_brick(brick)
    
    for boundary_k in range(len(boundaries)):
        for control_type_k in range(len(control_types)):
            name = "Control "+control_types[control_type_k]+" "+boundaries[boundary_k]
            expression = control_expressions[boundary_k][control_type_k]
            inse.set_control(name, expression)
    
    inse.set_initial_value("psi", psi_0)
    
    inse.set_time_scheme(t_0=0., t_f=t_f, dt=dt, dt_save=dt_save, init_step=init_step,
                         ts_adapt_time_step_increase_delay=0,
                         ts_adapt_clip_low=0.9, ts_adapt_clip_high=1.1,
                         ts_adapt_dt_min=ts_adapt_dt_min, ts_adapt_dt_max=dt_save,
                         ts_type=ts_type, bdf_order=bdf_order, ksp_type=ksp_type,
                         pc_type="lu", pc_factor_mat_solver_type=solver_type
                         )
    
    inse.solve()
    
    inse.hamiltonian.add_term(S.Term("Kinetic energy", "rho * psi * e_d / mu", [1]))
    inse.plot_Hamiltonian(with_powers=False)
    
    # psi = inse.get_solution("psi")
    # e_d = inse.get_solution("e_d")
    
    inse.export_to_pv("psi")
    inse.export_to_pv("e_d")
    
    print("Number of dofs:", inse.gf_model.nbdof())
    
    return inse

if __name__ == '__main__':

    inse = INSE()