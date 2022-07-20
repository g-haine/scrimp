#!/usr/bin/env python
# coding: utf-8

# For more informations about **SCRIMP**, see our [website](https://g-haine.github.io/scrimp/).
# 
# # The lossless 2D wave equation with Neumann boundary control
# 
# This notebook is meant to present the 2D lossless wave equation as first example of the SCRIMP wrapper for PFEM.
# 
# ## The model
# 
# Let us consider the vertical deflection from equilibrium $w$ of a 2D membrane $\Omega \subset \mathbb{R}^2$. Denoting $\rho$ the mass density and $T$ the Young modulus of the membrane, a positive definite tensor, leads to the following well-known *wave equation*
# $$
# \rho(x) \frac{\partial^2}{\partial t^2} w(t,x) - {\rm div} \left( T(x) \cdot {\rm grad} \left( w(t,x) \right) \right) = 0, \quad t \ge 0, \, x \in \Omega,
# $$
# together with *Neumann boundary control* 
# $$
# \left( T(x) \cdot {\rm grad} \left( w(t,x) \right) \right) \cdot \mathbf{n} = u_\partial(t,x), \quad t \ge 0, \, x \in \partial \Omega,
# $$
# where $\mathbf{n}$ is the outward normal to $\Omega$.
# 
# The **Hamiltonian** is the total mechanical energy, given as the sum of potential and kinetic energies
# $$
# \mathcal{H}(t) := \frac{1}{2} \int_\Omega \left( {\rm grad} \left( w(t,x) \right) \right)^\top \cdot T(x) \cdot {\rm grad} \left( w(t,x) \right) {\rm d}x 
# + \frac{1}{2} \int_\Omega \rho(x) \left( \frac{\partial}{\partial t} w(t,x) \right)^2 {\rm d}x, \qquad t \ge 0.
# $$
# Taking the *strain* $\mathbf{\alpha}_q := {\rm grad} \left( w \right)$ and the *linear momentum* $\alpha_p := \frac{\partial}{\partial t} w$ as **energy variables**, the Hamiltonian rewrites
# $$
# \mathcal{H}(t) = \mathcal{H}(\mathbf{\alpha}_q(t,\cdot), \alpha_p(t,\cdot)) = \frac{1}{2} \int_\Omega \left( \mathbf{\alpha}_q(t,x) \right)^\top \cdot T(x) \cdot \mathbf{\alpha}_q(t,x) {\rm d}x
# + \frac{1}{2} \int_\Omega \frac{\alpha_p^2(t,x)}{\rho(x)} {\rm d}x.
# $$
# The **co-energy variables** are by definition the variational derivatives of the Hamiltonian
# $$
# \mathbf{e}_q := \delta_{\mathbf{\alpha}_q} \mathcal{H} = T \cdot \mathbf{\alpha}_q, 
# \qquad e_p := \delta_{\alpha_p} \mathcal{H} = \frac{\alpha_p}{\rho},
# $$
# *i.e.* the *stress* and the *velocity* respectively. These equality are the **constitutive relation** which close the system.
# 
# Thanks to these variables, the wave equation writes as a **port-Hamiltonian system**
# $$
# \begin{pmatrix}
# \frac{\partial}{\partial t} \mathbf{\alpha}_q \\
# \frac{\partial}{\partial t} \alpha_p
# \end{pmatrix}
# =
# \begin{bmatrix}
# 0 & {\rm grad} \\
# {\rm div} & 0
# \end{bmatrix}
# \begin{pmatrix}
# \mathbf{e}_q \\
# e_p
# \end{pmatrix}, 
# \qquad \left\lbrace
# \begin{array}{rcl}
# \mathbf{e}_q &=& T \cdot \mathbf{\alpha}_q, \\
# e_p &=& \frac{\alpha_p}{\rho},
# \end{array}\right.
# \qquad \left\lbrace
# \begin{array}{rcl}
# u_\partial &=& \mathbf{e}_q \cdot \mathbf{n}, \\
# y_\partial &=& e_p|_{\partial \Omega},
# \end{array}\right.
# $$
# 
# The **power balance** satisfied by the Hamiltonian is
# $$
# \frac{\rm d}{{\rm d}t} \mathcal{H} = \langle u_\partial, y_\partial \rangle_{H^{-\frac12}(\partial \Omega),H^\frac12(\partial \Omega)}
# $$
# 
# To get rid of the algebraic constraints induced by the constitutive relations, one rewrites the port-Hamiltonian system as
# $$
# \begin{bmatrix}
# T^{-1} & 0 \\
# 0 & \rho
# \end{bmatrix}
# \begin{pmatrix}
# \frac{\partial}{\partial t} \mathbf{e}_q \\
# \frac{\partial}{\partial t} e_p
# \end{pmatrix}
# =
# \begin{bmatrix}
# 0 & {\rm grad} \\
# {\rm div} & 0
# \end{bmatrix}
# \begin{pmatrix}
# \mathbf{e}_q \\
# e_p
# \end{pmatrix}, 
# \qquad \left\lbrace
# \begin{array}{rcl}
# u_\partial &=& \mathbf{e}_q \cdot \mathbf{n}, \\
# y_\partial &=& e_p|_{\partial \Omega},
# \end{array}\right.
# $$
# known as the **co-energy formulation**. This allows to get a simple Ordinary Differential Equation at the discrete level (instead of a Differential Algebraic Equation in general).
# 
# ## The Partitioned Finite Element Method
# 
# The strategy follows three steps, inspired by the Mixed Finite Element Method for steady-state problem with homogeneous boundary condition
# * write the weak form of the system;
# * integrate by parts a **partition** of the state (such that $u_\partial$ appears); and
# * project on finite element spaces.
# 
# ### Weak formulation
# 
# Let $\phi_q$, $\varphi_p$ and $\psi$ be vector-valued, scalar-valued and boundary scalar-valued test functions respectively. The weak formulation reads
# $$
# \left\lbrace
# \begin{array}{rcl}
# \displaystyle \int_\Omega \phi_q \cdot T^{-1} \cdot \frac{\partial}{\partial t} \mathbf{e}_q 
# &=& \displaystyle \int_\Omega \phi_q \cdot {\rm grad} \left( e_p \right), \\
# \displaystyle \int_\Omega \varphi_p \rho \frac{\partial}{\partial t} e_p 
# &=& \displaystyle \int_\Omega \varphi_p {\rm div} \left( \mathbf{e}_q \right), \\
# \displaystyle \int_{\partial \Omega} \psi y_\partial &=& \displaystyle \int_{\partial \Omega} \psi e_p.
# \end{array}\right.
# $$
# 
# ### Integration by parts
# 
# The integration by parts of the second line makes $u_\partial = \mathbf{e}_q \cdot \mathbf{n}$ appear
# $$
# \left\lbrace
# \begin{array}{rcl}
# \displaystyle \int_\Omega \phi_q \cdot T^{-1} \cdot \frac{\partial}{\partial t} \mathbf{e}_q 
# &=& \displaystyle \int_\Omega \phi_q \cdot {\rm grad} \left( e_p \right), \\
# \displaystyle \int_\Omega \varphi_p \rho \frac{\partial}{\partial t} e_p 
# &=& \displaystyle - \int_\Omega {\rm grad} \left( \varphi_p \right) \cdot \mathbf{e}_q
# + \int_{\partial \Omega} \varphi_p u_\partial, \\
# \displaystyle \int_{\partial \Omega} \psi y_\partial &=& \displaystyle \int_{\partial \Omega} \psi e_p.
# \end{array}\right.
# $$
# 
# ### Projection
# 
# Let $(\phi_q^i)_{1 \le i \le N_q}$, $(\varphi_p^j)_{1 \le j \le N_p}$ and $(\psi^k)_{1 \le k \le N_\partial}$ be finite element families for $q$-type, $p$-type and boundary-type variables. Variables are approximated in their respective finite element family
# $$
# \mathbf{e}_q^d(t,x) := \sum_{i=1}^{N_q} e_q^i(t) \phi_q^i(x),
# \qquad e_p^d(t,x) := \sum_{j=1}^{N_p} e_p^j(t) \varphi_p^j(x),
# $$
# $$
# u_\partial^d(t,x) := \sum_{k=1}^{N_\partial} u_\partial^k(t) \psi^k(x),
# \qquad y_\partial^d(t,x) := \sum_{k=1}^{N_\partial} y_\partial^k(t) \psi^k(x).
# $$
# Denoting $\underline{\star}$ the (time-varying) vector of coordinates of the discretisation $\star^d$ of $\star$ in its respective finite element family, the discrete system reads
# $$
# \begin{bmatrix}
# M_q & 0 & 0 \\
# 0 & M_p & 0 \\
# 0 & 0 & M_\partial
# \end{bmatrix}
# \begin{pmatrix}
# \frac{\rm d}{{\rm d}t} \underline{e_q}(t) \\
# \frac{\rm d}{{\rm d}t} \underline{e_p}(t) \\
# - \underline{y_\partial}(t)
# \end{pmatrix}
# =
# \begin{bmatrix}
# 0 & D & 0 \\
# -D^\top & 0 & B \\
# 0 & -B^\top & 0
# \end{bmatrix}
# \begin{pmatrix}
# \underline{e_q}(t) \\
# \underline{e_p}(t) \\
# \underline{u_\partial}(t)
# \end{pmatrix}
# $$
# where
# $$
# (M_q)_{ij} := \int_\Omega \phi_q^i \cdot T^{-1} \cdot \phi_q^j,
# \qquad 
# (M_p)_{ij} := \int_\Omega \varphi_p^i \rho \varphi_p^j,
# \qquad 
# (M_\partial)_{ij} := \int_{\partial \Omega} \psi^i \psi^j,
# $$
# and
# $$
# (D)_{ij} := \int_\Omega \phi_q^i \cdot {\rm grad} \left( \varphi_p^j \right),
# \qquad
# (B)_{jk} := \int_{\partial \Omega} \varphi_p^j \psi^k.
# $$
# 
# ### Discrete Hamiltonian
# 
# By definition, the discrete Hamiltonian is equal to the continuous Hamiltonian evaluated in the approximated variables. As we are working with the co-energy formulation, a first step is to restate the Hamiltonian in terms of co-energy variables
# $$
# \mathcal{H} = \frac{1}{2} \int_\Omega \mathbf{e}_q \cdot T^{-1} \cdot \mathbf{e}_q 
# + \frac{1}{2} \int_\Omega \rho (e_p)^2.
# $$
# Then, the discrete Hamiltonian is defined as
# $$
# \mathcal{H}^d := \frac{1}{2} \int_\Omega \mathbf{e}_q^d \cdot T^{-1} \cdot \mathbf{e}_q^d 
# + \frac{1}{2} \int_\Omega \rho (e_p^d)^2.
# $$
# After straightforward computations, it comes
# $$
# \mathcal{H}^d(t) = \frac{1}{2} \underline{e_q}(t)^\top M_q \underline{e_q}(t) + \frac{1}{2} \underline{e_p}(t)^\top M_p \underline{e_p}(t),
# $$
# and the **discrete power balance** follows
# $$
# \frac{\rm d}{{\rm d}t} \mathcal{H}^d(t) = \underline{u_\partial}(t)^\top M_\partial \underline{y_\partial}(t).
# $$
# 

# ## 1: Load SCRIMP

import os
# Get the location of the current notebook
path = os.getcwd().split(os.sep)
# Guess the path to SCRIMP (assuming this notebook has not been moved)
rootdir = ""
for dir in path[0:-3]:
    rootdir = rootdir+dir+os.sep
print("Location of root directory for SCRIMP is:", rootdir)

# Append SCRIMP root directory to path
import sys
sys.path.append(rootdir)

# Append all sub-directories to path
import scrimp
scrimp.setup_path(rootdir)


# ## 2: Load third-party libraries

import numpy as np # Array manipulations
import matplotlib.pyplot as plt # To plot figures
import fenics as fe # Finite elements library
import multiphenics as mpfe # Allow easier use of FEniCS for Mixed Finite Elements spaces

# Can be found in src/modules directory of SCRIMP
import gmsh_utils
import meshio_utils
import fenics_utils
import multiphenics_utils

fe.SubSystemsManager.init_petsc()
from petsc4py import PETSc
import PETSc_utils


# ## 3: Define geometry from ``.geo`` file and mesh with GMSH

# May be any 2D .geo file in data/geo_file
Geometry = "Rectangle"
# The following parameters are those from .geo files (see data/geo_files directory to adapt them)
L = 1.0
l = 0.5
h = 0.05
params = {'h':h,'hmin':h,'layer':1./2.5,'L':L,'l':l}

# Define string paths
geofile = rootdir+"data/geo_files/"+Geometry+".geo"
gmshfile = rootdir+"data/geo_files/"+Geometry+".msh"

# Dimension must be given explicitly
dim = 2

# Mesh generation
# /!\ : Dolfin (FEniCS) only fully supports first-order meshes
gmsh_utils.generate_mesh_cli(geofile,gmshfile,dim,refinement=1,log=1,                             parameters=params,order=1,gmshformat=2,binary=True)
# Print mesh informations
gmsh_utils.print_summary(gmshfile)


# ## 4: Load XDMF mesh for FEniCS

# Check .msh file
meshio_utils.print_diagnostic_gmsh(gmshfile)
# Convert to XDMF file
xdmfiles = meshio_utils.convert_gmsh_to_XDMF(gmshfile,prune_z=True) # prune_z prunes the third coordinate
# Load in Dolfin format (FEniCS)
dmesh = fenics_utils.DolfinMesh.init_from_xdmf(xdmfiles['tetra'*(dim==3)+'triangle'*(dim==2)],
                                               xdmfiles['triangle'*(dim==3)+'line'*(dim==2)])
fe.plot(dmesh.mesh)

# Get the physical tags from .geo file
meshtag = meshio_utils.read_phys_gmsh(xdmfiles['gmsh_physical_entities'])
# Boundaries
bnd_name = "Gamma" # From .geo file, may be a list
Gamma = meshtag[dim-2][bnd_name][0]
# Restriction of the mesh at the boundaries, useful for the space definition
Gamma_restriction = multiphenics_utils.build_MeshRestriction_from_tags(dmesh.mesh,
                                                                       dmesh.boundaries,[Gamma])


# ## 5: Finite dimensional spaces

# FE for \alpha_q and e_q
type_q = "RT"
order_q = 2
# FE for \alpha_p and e_p
type_p = "CG"
order_p = 3
# FE for u and y
type_b = "CG"
order_b = 3

# Care must be taken: for Galerkin FE, VectorFunctionSpace must be used instead of FunctionSpace
if ((type_q == 'CG') or (type_q == 'DG')):
    V_q = fe.VectorFunctionSpace(dmesh.mesh, type_q, order_q)
else:
    V_q = fe.FunctionSpace(dmesh.mesh, type_q, order_q)
V_p = fe.FunctionSpace(dmesh.mesh, type_p, order_p)
V_u = fe.FunctionSpace(dmesh.mesh, type_b, order_b)
V_y = fe.FunctionSpace(dmesh.mesh, type_b, order_b)

# Block function space (Mixed Finite Elements)
W = mpfe.BlockFunctionSpace([V_q, V_p, V_u, V_y], 
                            restrict=[None, None, Gamma_restriction, Gamma_restriction])

# Test/Trial functions
fun_test = mpfe.block_split(mpfe.BlockTestFunction(W))
fun_trial = mpfe.block_split(mpfe.BlockTrialFunction(W))
                                


# ## 6: Dirac structure and Jacobian

# Young modulus and its inverse (for the sake of simplicity)
Txx = "2."
Tyy = "1."
Txy = "0.2 * (1+x[0]) * (1-x[0])"
T = fe.Expression([[Txx,Txy],
                   [Txy,Tyy]], degree=2)
Tinv = fe.Expression([[Tyy,"-"+Txy],
                      ["-"+Txy,Txx]], degree=2)/fe.det(T)

# Mass density
rho = fe.Expression("2. + 0.25 * (1+x[0]) * (1-x[0])", degree=2)

def Construct_Matrices(fun_trial,fun_test,dx,ds):
    
    # Split trial/test functions
    (e_q, e_p, u_b, y_b) = fun_trial
    (v_q, v_p, v_u, v_y) = fun_test
    
    # Mass matrices
    M_q = fe.dot(e_q,fe.dot(Tinv,v_q))*dx
    M_p = e_p*rho*v_p*dx
    M_u = u_b*v_u*ds
    M_y = y_b*v_y*ds
    
    # Matrices constituting J
    D = fe.dot(fe.grad(e_p), v_q)*dx
    Dt = fe.dot(e_q, fe.grad(v_p))*dx
    B = u_b*v_p*ds
    Bt = e_p*v_y*ds
    
    return M_q, M_p, M_u, M_y, D, Dt, B, Bt

M_q, M_p, M_u, M_y, D, Dt, B, Bt = Construct_Matrices(fun_trial,fun_test,dmesh.dx,dmesh.ds)

def Dirac(M_q, M_p, M_u, M_y, D, Dt, B, Bt):
    
    # It is not necessary to add null variational formulations, multiphenics handles this
    # Matrices M and J are rewritten with E and A as : E dz/dt = A z
    E = [
            # - M_q d(e_q)/dt
            [ -M_q ],
            # - M_p d(e_p)/dt
            [ -M_p ]
        ]
    
    A = [
            # D e_p
            [ D ],
            # -D^T e_q + B u_b
            [ -Dt, B ],
            # M_y y = -B^T e_p
            [ -M_y, -Bt ],
        
            # The following will allow to treat time-varying control "outside" the Dirac method
            # using a method defining only u_b thanks to v_b: see Define_Control method
            [ -M_u ]
        ]
    
    return (E, A)

# Jacobian consists of E and A when considering the general DAE form: F(t,z,dz/dt) = 0
Jacobian = lambda : Dirac(M_q, M_p, M_u, M_y, D, Dt, B, Bt)


# ## 7: Boundary source term

# UFL Format
t_control_begin = 0.5
t_control_end = 2.5
control = "5. * x[0] * x[0] * x[1] * sin(t-"+str(t_control_begin)+") * (t >"+str(t_control_begin)             +") * sin("+str(t_control_end)+"-t) * (t < "+str(t_control_end)+")"
u = fe.Expression(control, element=V_u.ufl_element(), domain=dmesh.mesh, t=0.)

def Define_Control(t,u,fun_test,ds):
    
    u.t = t # u at time t
    
    (v_q, v_p, v_u, v_y) = fun_test # Split test functions
    
    return [ [ u*v_u*ds ] ] # with the last line of the Dirac method, this defines u_b with u

# Residual consists of the time-varying forcing when considering the general DAE form: F(t,z,dz/dt) = 0
Residual = lambda t : Define_Control(t,u,fun_test,dmesh.ds)


# # 8: Initial data

x0Array = np.zeros(W.dim(),) # Initialisation q0 and p0

q0 = fe.Expression(["0.","0."], element=V_q.ufl_element(), domain=dmesh.mesh)
q0Fun = fe.interpolate(q0, V_q) # Warning: compatible condition "initial data -- boundary control" must be fulfilled
x0Array[0:V_q.dim()] = np.array(q0Fun.vector())

p0 = fe.Expression("2.*exp(-50.*((x[0]-L/2.)*(x[0]-L/2.)+(x[1]-l/2.)*(x[1]-l/2.)))", element=V_p.ufl_element(), domain=dmesh.mesh,L=L,l=l)
p0Fun = fe.interpolate(p0, V_p)
x0Array[V_q.dim():V_q.dim()+V_p.dim()] = np.array(p0Fun.vector())

x0 = PETSc.Vec().createWithArray(x0Array)


# # 9: Time resolution

# PETSc TS solver options (time integration)
# OptDB collects options for PETSc TS: https://petsc.org/release/docs/manualpages/TS/index.html
# a statement as -option can be passed to PETSc TS by setting OptDB["option"]
OptDB = PETSc_utils.get_cleared_options_db()
PETSc.Sys.pushErrorHandler("python") # To see correctly the errors: interferences between FEniCS and PETSc???

OptDB["ts_type"] = "cn" # Crank-Nicolson
OptDB["ksp_type"] = "preonly" # Direct linear solver
OptDB["pc_type"] = "lu" # LU preconditioner

dt = 0.01
tf = 5.
Integration_params = {'dt_export':5.*dt}

# u_b and y_b are algebraic
idx_Lagrange_multipliers = np.concatenate((W.block_dofmap().block_owned_dofs__global_numbering(2),
                                           W.block_dofmap().block_owned_dofs__global_numbering(3)))

name = "2D Wave equation -- Neumann boundary control"
dae = multiphenics_utils.PDAE_linear(Residual,Jacobian,
                                     W.dim(),name,idx_alg=idx_Lagrange_multipliers)

PETSc_utils.TS_integration_dae(dae,x0,dt,tf,**Integration_params)

t_sol = dae.history['t']
z_sol = dae.history['z']


# # 10: Post-processing

def Plot_Hamiltonian(t,z,u,Tinv,rho,fun_trial,fun_test,dx,ds):
    
    # Number of time steps
    Nt = np.array(t).size
    
    # For the potential energy
    Hamq = np.zeros((Nt,))
    e_q = mpfe.BlockFunction(W,mpfe.la.BlockPETScVector(z[0]))[0]
    Hamq[0] = 0.5*fe.assemble(fe.dot(e_q,fe.dot(Tinv, e_q))*dx)
    
    # For the kinetic energy
    Hamp = np.zeros((Nt,))
    e_p = mpfe.BlockFunction(W,mpfe.la.BlockPETScVector(z[0]))[1]
    Hamp[0] = 0.5*fe.assemble(e_p*rho*e_p*dx)
    
    # For the total mechanical energy (the Hamiltonian)
    Ham = np.zeros((Nt,))
    Ham[0] = Hamq[0] + Hamp[0]
    
    # For supplied energy
    Supplied = np.zeros((Nt,))
    
    for k in range(1,Nt):
        # Extract solutions
        e_q = mpfe.BlockFunction(W,mpfe.la.BlockPETScVector(z[k]))[0]
        e_p = mpfe.BlockFunction(W,mpfe.la.BlockPETScVector(z[k]))[1]
        y_b_old = mpfe.BlockFunction(W,mpfe.la.BlockPETScVector(z[k-1]))[3]
        y_b = mpfe.BlockFunction(W,mpfe.la.BlockPETScVector(z[k]))[3]
        # Compute energies
        Hamq[k] = 0.5*fe.assemble(fe.dot(e_q,fe.dot(Tinv, e_q))*dx)
        Hamp[k] = 0.5*fe.assemble(e_p*rho*e_p*dx)
        Ham[k] = Hamq[k] + Hamp[k]
        u.t = t[k-1]
        Supplied[k] = Supplied[k-1] + 0.5 * (t[k] - t[k-1]) * fe.assemble(y_b_old*u*ds)
        u.t = t[k]
        Supplied[k] = Supplied[k] + 0.5 * (t[k] - t[k-1]) * fe.assemble(y_b*u*ds)
    
    fig = plt.figure(figsize=(16,9))
    plt.plot(t,Ham+Supplied,t,Ham,t,Hamq,t,Hamp,t,-Supplied)
    plt.title("Wave 2D -- Neuman boundary control ("
              +str(W.tabulate_dof_coordinates().shape[0])+" dofs)")
    plt.xlabel("Time [s]")
    plt.ylabel("Energies [J]")
    plt.grid(axis="both")
    plt.legend(["Total energy","Hamiltonian","Potential","Kinetic","Supplied"])
    plt.show()

Plot_Hamiltonian(t_sol,z_sol,u,Tinv,rho,fun_trial,fun_test,dmesh.dx,dmesh.ds)


# Plot some instants (in % of final time)
instants = (0.,0.25,0.5,0.75,1.)

# Get coordinates
dofs = V_p.tabulate_dof_coordinates().reshape((-1, dmesh.mesh.geometry().dim()))
x = [dof[0] for dof in dofs]
y = [dof[1] for dof in dofs]

def Plot_Instant(t,z,W,V_p,instant):
    idx = min(int(instant*len(t)),len(t)-1)
    e_q = mpfe.BlockFunction(W,mpfe.la.BlockPETScVector(z[idx]))[0]
    e_p = mpfe.BlockFunction(W,mpfe.la.BlockPETScVector(z[idx]))[1]
    
    fig = plt.figure(figsize=(16,9))
    ax1 = fig.add_subplot(2, 2, 1)
    plt.title("Wave equation -- Stress at time t="+str(t[idx])+" s")
    fe.plot(e_q)
    ax2 = fig.add_subplot(2, 2, 2, projection='3d')
    plt.title("Wave equation -- Velocity at time t="+str(t[idx])+" s")
    ax2.plot_trisurf(x,y,e_p.vector(), cmap=plt.cm.coolwarm)
    ax2.set_zlim(-10, 10)
    plt.show()

for instant in instants:
    Plot_Instant(t_sol,z_sol,W,V_p,instant)


# # 11: Exports

# Export for ParaView
pvdfile = 'vtk/'
sol_idx=[0,1]; sol_name = ["e_q","e_p"]
print(f"Export for paraview...")
multiphenics_utils.export_pvd(pvdfile,sol_idx,sol_name,t_sol,z_sol,W)
print("Done.")

# Export matrices to MATLAB (binary) format
from scipy.io import savemat
print(f"Export for MATLAB...")
Mdic = {"M_q": fe.assemble(M_q), "M_p": fe.assemble(M_p), 
        "M_u": fe.assemble(M_u), "D": fe.assemble(D), 
        "B": fe.assemble(B)}
savemat("wave_matrices.mat", Mdic)
print("Done.")


# Export this .ipynb to a simple .py script
get_ipython().system('jupyter nbconvert --to script --no-prompt wave.ipynb')


# Export this .ipynb to markdown
get_ipython().system('jupyter nbconvert --to markdown --no-prompt wave.ipynb')

