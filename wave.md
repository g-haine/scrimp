# The 2D wave equation with Neumann boundary control

Let us consider the vertical deflection from equilibrium $w$ of a 2D membrane $\Omega \subset \mathbb{R}^2$. Denoting $\rho$ the mass density and $T$ the Young modulus of the membrane, a positive definite tensor, leads to the following well-known *wave equation*

$$
\rho(x) \frac{\partial^2}{\partial t^2} w(t,x) - {\rm div} \left( T(x) \cdot {\rm grad} \left( w(t,x) \right) \right) = 0, \quad t \ge 0, \, x \in \Omega,
$$

together with *Neumann boundary control* 
$$
\left( T(x) \cdot {\rm grad} \left( w(t,x) \right) \right) \cdot \mathbf{n} = u_\partial(t,x), \quad t \ge 0, \, x \in \partial \Omega,
$$
where $\mathbf{n}$ is the outward normal to $\Omega$.

The **Hamiltonian** is the total mechanical energy, given as the sum of potential and kinetic energies
$$
\mathcal{H}(t) := \frac{1}{2} \int_\Omega \left( {\rm grad} \left( w(t,x) \right) \right)^\top \cdot T(x) \cdot {\rm grad} \left( w(t,x) \right) {\rm d}x 
+ \frac{1}{2} \int_\Omega \rho(x) \left( \frac{\partial}{\partial t} w(t,x) \right)^2 {\rm d}x, \qquad t \ge 0.
$$
Taking the *strain* $\mathbf{\alpha}_q := {\rm grad} \left( w \right)$ and the *linear momentum* $\alpha_p := \frac{\partial}{\partial t} w$ as **energy variables**, the Hamiltonian rewrites
$$
\mathcal{H}(t) = \mathcal{H}(\mathbf{\alpha}_q(t,\cdot), \alpha_p(t,\cdot)) = \frac{1}{2} \int_\Omega \left( \mathbf{\alpha}_q(t,x) \right)^\top \cdot T(x) \cdot \mathbf{\alpha}_q(t,x) {\rm d}x
+ \frac{1}{2} \int_\Omega \frac{\alpha_p^2(t,x)}{\rho(x)} {\rm d}x.
$$
The **co-energy variables** are by definition the variational derivatives of the Hamiltonian
$$
\mathbf{e}_q := \delta_{\mathbf{\alpha}_q} \mathcal{H} = T \cdot \mathbf{\alpha}_q, 
\qquad e_p := \delta_{\alpha_p} \mathcal{H} = \frac{\alpha_p}{\rho},
$$
*i.e.* the *stress* and the *velocity* respectively. These equality are the **constitutive relation** which close the system.

Thanks to these variables, the wave equation writes as a **port-Hamiltonian system**
$$
\begin{pmatrix}
\frac{\partial}{\partial t} \mathbf{\alpha}_q \\
\frac{\partial}{\partial t} \alpha_p
\end{pmatrix}
=
\begin{bmatrix}
0 & {\rm grad} \\
{\rm div} & 0
\end{bmatrix}
\begin{pmatrix}
\mathbf{e}_q \\
e_p
\end{pmatrix}, 
\qquad \left\lbrace
\begin{array}{rcl}
\mathbf{e}_q &=& T \cdot \mathbf{\alpha}_q, \\
e_p &=& \frac{\alpha_p}{\rho},
\end{array}\right.
\qquad \left\lbrace
\begin{array}{rcl}
u_\partial &=& \mathbf{e}_q \cdot \mathbf{n}, \\
y_\partial &=& e_p|_{\partial \Omega},
\end{array}\right.
$$

The **power balance** satisfied by the Hamiltonian is
$$
\frac{\rm d}{{\rm d}t} \mathcal{H} = \langle u_\partial, y_\partial \rangle_{H^{-\frac12}(\partial \Omega),H^\frac12(\partial \Omega)}
$$

To get rid of the algebraic constraints induced by the constitutive relations, one rewrites the port-Hamiltonian system as
$$
\begin{bmatrix}
T^{-1} & 0 \\
0 & \rho
\end{bmatrix}
\begin{pmatrix}
\frac{\partial}{\partial t} \mathbf{e}_q \\
\frac{\partial}{\partial t} e_p
\end{pmatrix}
=
\begin{bmatrix}
0 & {\rm grad} \\
{\rm div} & 0
\end{bmatrix}
\begin{pmatrix}
\mathbf{e}_q \\
e_p
\end{pmatrix}, 
\qquad \left\lbrace
\begin{array}{rcl}
u_\partial &=& \mathbf{e}_q \cdot \mathbf{n}, \\
y_\partial &=& e_p|_{\partial \Omega},
\end{array}\right.
$$
known as the **co-energy formulation**. This allows to get a simple Ordinary Differential Equation at the discrete level (instead of a Differential Algebraic Equation in general).

# The Partitioned Finite Element Method

The strategy follows three steps, inspired by the Mixed Finite Element Method for steady-state problem with homogeneous boundary condition
* write the weak form of the system;
* integrate by parts a **partition** of the state (such that $u_\partial$ appears); and
* project on finite element spaces.

## Weak formulation

Let $\phi_q$, $\varphi_p$ and $\psi$ be vector-valued, scalar-valued and boundary scalar-valued test functions respectively. The weak formulation reads
$$
\left\lbrace
\begin{array}{rcl}
\displaystyle \int_\Omega \phi_q \cdot T^{-1} \cdot \frac{\partial}{\partial t} \mathbf{e}_q 
&=& \displaystyle \int_\Omega \phi_q \cdot {\rm grad} \left( e_p \right), \\
\displaystyle \int_\Omega \varphi_p \rho \frac{\partial}{\partial t} e_p 
&=& \displaystyle \int_\Omega \varphi_p {\rm div} \left( \mathbf{e}_q \right), \\
\displaystyle \int_{\partial \Omega} \psi y_\partial &=& \displaystyle \int_{\partial \Omega} \psi e_p.
\end{array}\right.
$$

## Integration by parts

The integration by parts of the second line makes $u_\partial = \mathbf{e}_q \cdot \mathbf{n}$ appear
$$
\left\lbrace
\begin{array}{rcl}
\displaystyle \int_\Omega \phi_q \cdot T^{-1} \cdot \frac{\partial}{\partial t} \mathbf{e}_q 
&=& \displaystyle \int_\Omega \phi_q \cdot {\rm grad} \left( e_p \right), \\
\displaystyle \int_\Omega \varphi_p \rho \frac{\partial}{\partial t} e_p 
&=& \displaystyle - \int_\Omega {\rm grad} \left( \varphi_p \right) \cdot \mathbf{e}_q
+ \int_{\partial \Omega} \varphi_p u_\partial, \\
\displaystyle \int_{\partial \Omega} \psi y_\partial &=& \displaystyle \int_{\partial \Omega} \psi e_p.
\end{array}\right.
$$

## Projection

Let $(\phi_q^i)_{1 \le i \le N_q}$, $(\varphi_p^j)_{1 \le j \le N_p}$ and $(\psi^k)_{1 \le k \le N_\partial}$ be finite element families for $q$-type, $p$-type and boundary-type variables. Variables are approximated in their respective finite element family
$$
\mathbf{e}_q^d(t,x) := \sum_{i=1}^{N_q} e_q^i(t) \phi_q^i(x),
\qquad e_p^d(t,x) := \sum_{j=1}^{N_p} e_p^j(t) \varphi_p^j(x),
$$
$$
u_\partial^d(t,x) := \sum_{k=1}^{N_\partial} u_\partial^k(t) \psi^k(x),
\qquad y_\partial^d(t,x) := \sum_{k=1}^{N_\partial} y_\partial^k(t) \psi^k(x).
$$
Denoting $\underline{\star}$ the (time-varying) vector of coordinates of the discretisation $\star^d$ of $\star$ in its respective finite element family, the discrete system reads
$$
\begin{bmatrix}
M_q & 0 & 0 \\
0 & M_p & 0 \\
0 & 0 & M_\partial
\end{bmatrix}
\begin{pmatrix}
\frac{\rm d}{{\rm d}t} \underline{e_q}(t) \\
\frac{\rm d}{{\rm d}t} \underline{e_p}(t) \\
- \underline{y_\partial}(t)
\end{pmatrix}
=
\begin{bmatrix}
0 & D & 0 \\
-D^\top & 0 & B \\
0 & -B^\top & 0
\end{bmatrix}
\begin{pmatrix}
\underline{e_q}(t) \\
\underline{e_p}(t) \\
\underline{u_\partial}(t)
\end{pmatrix}
$$
where
$$
(M_q)_{ij} := \int_\Omega \phi_q^i \cdot T^{-1} \cdot \phi_q^j,
\qquad 
(M_p)_{ij} := \int_\Omega \varphi_p^i \rho \varphi_p^j,
\qquad 
(M_\partial)_{ij} := \int_{\partial \Omega} \psi^i \psi^j,
$$
and
$$
(D)_{ij} := \int_\Omega \phi_q^i \cdot {\rm grad} \left( \varphi_p^j \right),
\qquad
(B)_{jk} := \int_{\partial \Omega} \varphi_p^j \psi^k.
$$

## Discrete Hamiltonian

By definition, the discrete Hamiltonian is equal to the continuous Hamiltonian evaluated in the approximated variables. As we are working with the co-energy formulation, a first step is to restate the Hamiltonian in terms of co-energy variables
$$
\mathcal{H} = \frac{1}{2} \int_\Omega \mathbf{e}_q \cdot T^{-1} \cdot \mathbf{e}_q 
+ \frac{1}{2} \int_\Omega \rho (e_p)^2.
$$
Then, the discrete Hamiltonian is defined as
$$
\mathcal{H}^d := \frac{1}{2} \int_\Omega \mathbf{e}_q^d \cdot T^{-1} \cdot \mathbf{e}_q^d 
+ \frac{1}{2} \int_\Omega \rho (e_p^d)^2.
$$
After straightforward computations, it comes
$$
\mathcal{H}^d(t) = \frac{1}{2} \underline{e_q}(t)^\top M_q \underline{e_q}(t) + \frac{1}{2} \underline{e_p}(t)^\top M_p \underline{e_p}(t),
$$
and the **discrete power balance** follows
$$
\frac{\rm d}{{\rm d}t} \mathcal{H}^d(t) = \underline{u_\partial}(t)^\top M_\partial \underline{y_\partial}(t).
$$
