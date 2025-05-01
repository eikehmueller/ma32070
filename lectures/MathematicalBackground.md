# Mathematical background
## Model problem
In this course we will focus on the following PDE of the diffusion-reaction type in some bounded domain $\Omega\subset \mathbb{R}^2$:
$$
-\nabla \cdot (\kappa \nabla  u(x)) + \omega\; u(x) = f(x) \qquad \text{for $x\in \Omega$}\qquad(\dagger)
$$
with boundary condition $\kappa\; n\cdot \nabla u(x)=g(x)$ for $x\in\partial \Omega$. We assume that $\omega, \kappa>0$ are positive constants and $f(x)$, $g(x)$ are given functions. Using zero-based indexing (as is used in Python) we will write $x=(x_0,x_1)\in\mathbb{R}^2$ such that $\nabla=(\frac{\partial}{\partial x_0},\frac{\partial}{\partial x_1})^\top$ is the nabla-operator. Note that in the case $\kappa=1$, $\omega=0$ the problem would reduce to the Poisson equation $-\Delta u(x)=f(x)$. Unfortunately, for the given boundary condition the solution of the Poisson equation is not unique, which is why we do not consider this case here. However, the methods developed in this course can be readily applied to this setup, provided we extend them to treat Dirichlet boundary conditions of the form $u(x)=\widetilde{g}(x)$ for $x\in\partial \Omega$.

## Weak solutions
To solve $(\dagger)$, we seek solutions $u(x)$ in some function space $V$. In the following we choose $V:=H^1(\Omega)\subset L_2(\Omega)$, which is the space of all real-valued functions on $\Omega$ which have a square-integrable first derivative. More specifically, define the following two norms
$$
\begin{aligned}
\| u\|_{L_2(\Omega)} &:= \left(\int_\Omega u(x)^2\;dx\right)^{\frac{1}{2}}\\
\| u\|_{V} = \| u\|_{H^1(\Omega)} &:= \left(\int_\Omega \left(u(x)^2+|\nabla u|^2\right)\;dx\right)^{\frac{1}{2}}
\end{aligned}
$$
and then set $L_2(\Omega) = \left\{u(x) : ||u||_{L_2(\Omega)}<\infty\right\}$ (the space of square-integrable real unctions) and $H^1(\Omega) = \left\{u(x) : ||u||_{H_1(\Omega)}<\infty\right\}$. Since in $(\dagger)$ two derivatives act on $u(x)$, we can only determine the solution in the **weak sense**: Find $u(x)\in V$ such that
$$
\int_\Omega \left(-v(x)\nabla \cdot(\kappa \nabla  u(x)) + \omega\; v(x) u(x)\right)\;dx = \int_\Omega f(x) v(x)\;dx \qquad \text{for all $v(x)\in V$}.
$$
Note that in contrast to $(\dagger)$ we no longer require that the equation is satisfied at every point $x$. Discussing in which sense these weak solutions are equivalent to solutions of $(\dagger)$ (which is sometimes also referred to as the **"strong"** form of the equation) is beyond the scope of this course. After integrating the first term under the integral on the left-hand side by parts, the weak form becomes
$$
\int_\Omega \left(\kappa \nabla v(x) \cdot \nabla  u(x) + \omega\; v(x) u(x)\right)\;dx - \int_{\partial \Omega } \kappa\;v(x) n\cdot \nabla(u)\;ds = \int_\Omega f(x) v(x)\;dx.
$$
Crucially, only first derivatives of the solution $u(x)$ and test function $v(x)$ are required now. Using the boundary condition $\kappa\; n\cdot \nabla u(x)=g(x)$ for $x\in\partial\Omega$, we can rewrite this as
$$
\int_\Omega \left(\kappa \nabla v(x) \cdot \nabla  u(x) + \omega\; v(x) u(x)\right)\;dx  = \int_\Omega f(x) v(x)\;dx + \int_{\partial \Omega} g(x) v(x)\;ds.
$$
Let us define the symmetric *bilinear form* $a(\cdot,\cdot): V\times V \rightarrow \mathbb{R}$  with
$$
a(u,v) := \int_\Omega \left(\kappa \nabla v(x) \cdot \nabla  u(x) + \omega\; v(x) u(x)\right)\;dx
$$
and the linear form $b(\cdot):V\rightarrow \mathbb{R}$ with
$$
b(v) := \int_\Omega f(x) v(x)\;dx+ \int_{\partial \Omega} g(x) v(x)\;ds.
$$
#### Exercise
Convince yourself that $a(\cdot,\cdot)$ and $b(\cdot)$ are indeed (bi-) linear:
* $a(c_1 u^{(1)} + c_2 u^{(2)}_2,v) = c_1 a(u^{(1)},v) + c_2 a(u^{(2)},v)$ for all $c_1,c_2\in \mathbb{R}$, $u^{(1)}, u^{(2)},v \in V$
* $a(u,c_1 v^{(1)} + c_2 v^{(2)}_2) = c_1 a(u,v^{(1)}) + c_2 a(u,v^{(2)})$ for all $c_1,c_2\in \mathbb{R}$, $u,v^{(1)}, v^{(2)} \in V$
* $b(c_1 v^{(1)} + c_2 v^{(2)}_2)=c_1b( v^{(1)}) + c_2 b(v^{(2)}_2)$ for all $c_1,c_2\in \mathbb{R}$, $v^{(1)}, v^{(2)} \in V$
  
and that $a(\cdot,\cdot)$ is symmetric:
* $a(v,u) = a(u,v)$ for all $u,v\in V$
  
With these (bi-)linear forms, we can formulate the weak problem as follows: Find $u(x)\in V$ such that
$$
a(u,v) = b(v) \qquad \text{for all $v(x)\in V$}.\qquad(\ddagger)
$$
## Finite element solutions
Now, obviously it is not possible to solve $(\ddagger)$ on a computer since $V$ contains infinitely many functions. Instead, we try to find solutions in a finite-dimensional subspace $V_h\subset V$. This could for example be the space of all functions that are piecewise linear on a given mesh with spacing $h$. We will be more precise about what that means later in this course. Then the problem becomes: find $u_h\in V_h$ such that 
$$
a(u_h,v_h) = b(v_h) \qquad \text{for all $v_h(x)\in V_h$ }.\qquad(\ddagger_h)
$$
### Existence and convergence of the solution
It can be shown that $(\ddagger)$ and $(\ddagger_h)$ have unique solutions provided the linear form $b(\cdot)$ and the bilinear form $a(\cdot,\cdot)$ satisfy the following two conditions:
* **Boundedness**: there exists some positive constant $C_+ > 0$ such that 
$$a(u,v) \le C_+ \|u\|_V \|v\|_V \qquad\text{and}$$
$$b(v) \le C_+ \|v\|_V \qquad\text{for all $u,v\in V$}.$$
* **Coercivity**: there exists some positive constant $C_- > 0$ such that
$$ 
a(u,u) \ge C_- \|u\|_V^2 \qquad\text{for all $u\in V$}.
$$
It turns out that both conditions are satisfied for the $a(\cdot,\cdot)$, $b(\cdot)$ defined above. Furthermore, the solutions satisfy $\|u\|_V,\|u_h\|_V\le C:=C_+/C$ and the difference between the solution $u_h(x)$ of $(\ddagger_h)$ and the solution $u(x)$ of $(\ddagger)$ can be bounded as follows:
$$
\|u_h - u\|_V \le C \min_{v_h\in V_h}\|u-v_h\|_V.
$$
The term $\min_{v_h\in V_h}\|u-v_h\|_V$ only depends on the choice of function spaces $V$, $V_h$ and describes how well the function $u(x) \in V$ can be approximated by a function $v_h\in V_h$. For a suitable choice of $V_h$, which we will discuss later, one can show that $\min_{v_h\in V_h}\|u-v_h\|_V\le C' h^{2\mu}$ for some positive integer $\mu\ge 1$ and positive constant $C'>0$. Hence, the finite element solution $u_h(x)$ converges to the "true'' solution $u(x)$ as the mesh is refined ($h\rightarrow 0$):
$$
\|u_h - u\|_V \le C C' h^{2\mu}.
$$
## Reduction to linear algebra problem
Since $V_h$ is finite dimensional, we can choose a basis $\{\Phi^{(h)}_j(x)\}_{j=0}^{n-1}$ such that every function $u_h(x)\in V_h$ can be written as
$$
u_h(x) = \sum_{j=0}^{n-1} u^{(h)}_j \Phi^{(h)}_j(x) \qquad\text{for all $x\in\Omega$}\qquad(\star).
$$
The vector $\boldsymbol{u}^{(h)}=(u^{(h)}_0,u^{(h)}_1,\dots,u^{(h)}_{n-1})\in\mathbb{R}^n$ is often referred to as the degrees-of-freedom vector (short: dof-vector) since its knowledge determines $u_h(x)$. Picking $v_h(x)=\Phi^{(h)}_i(x)$ and inserting the expansion of $u_h(x)$ in $(\star)$ into $(\ddagger_h)$ we obtain
$$
b^{(h)}_i:=b(\Phi^{(h)}_i) = a\left(\sum_{j=0}^{n-1} u^{(h)}_j \Phi^{(h)}_j,\Phi^{(h)}_i\right) = 
\sum_{j=0}^{n-1} u^{(h)}_j a\left( \Phi^{(h)}_i,\Phi^{(h)}_j\right),
$$
where we used the symmetry and bi-linearity of $a(\cdot,\cdot)$. Defining the vector $\boldsymbol{b}^{(h)} := (b(\Phi^{(h)}_0),b(\Phi^{(h)}_1,\dots,b(\Phi^{(h)}_{n-1})))$ and the $n\times n$ matrix $A^{(h)}$ with $A^{(h)}_{ij}:= a\left(\Phi^{(h)}_i,\Phi^{(h)}_j\right)$ we arrive at the following linear system for the dof-vector $\boldsymbol{u}^{(h)}$:
$$
A^{(h)} \boldsymbol{u}^{(h)} = \boldsymbol{b}^{(h)}.
$$
At this point it is worth stressing that although $\boldsymbol{u}^{(h)}$ and $\boldsymbol{b}^{(h)}$ are both vectors in $\mathbb{R}^n$, they are constructed in a fundamentally different way:

* The dof-vector $\boldsymbol{u}^{(h)}$ is a so-called **primal** vector: its components $U_j^{(h)}$ are the expansion coefficients of the function $u_h(x)$ in $(\star)$.
* In contrast, the right-hand-side vector $\boldsymbol{b}^{(h)}$ is a so-called **dual** vector: its components $b(\Phi_j^{(h)})$ are obtained by evaluating the linear functional $b(\cdot)$ for the basis functions.

The reason for this is that $b(\cdot)$ is an element of the dual space $V^*$, which consists of all linear functionals defined on the space $V$.
### Solution procedure
In summary, the solution procedure for $(\ddagger_h)$ is this:
1. Assemble the matrix $A^{(h)}$.
2. Assemble the right-hand-side vector $\boldsymbol{b}^{(h)}$.
3. Solve the linear system $A^{(h)} \boldsymbol{u}^{(h)} = \boldsymbol{b}^{(h)}$ for $\boldsymbol{u}^{(h)}$.
4. Reconstruct the solution $u_h(x)$ from the dof-vector $\boldsymbol{u}^{(h)}$ according to the expansion in $(\star)$.
