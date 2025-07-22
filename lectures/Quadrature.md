# Numerical quadrature
The weak form of the PDE is defined via suitable integrals such as $\int_\Omega v(x)f(x)\;dx$. In general, it is not possible to evaluate these integrals exactly. Furthermore, since the finite element discretisation (replacing $V\mapsto V_h$ and solving the associated linear algebra problem) already introduces an error, exact integration is not necessary, provided we can find an approximate integration method with errors that are of the same order of magnitude.

*&#169; Eike Mueller, University of Bath 2025. These lecture notes are copyright of Eike Mueller, University of Bath. They are provided exclusively for educational purposes at the University and are to be downloaded or copied for your private study only. Further distribution, e.g. by upload to external repositories, is prohibited.*


## Gauss-Legendre quadrature in one dimension
Numerical quadrature aims to approximate the integral of a function with a finite sum:

$$
\int_{-1}^{+1} f(z)\;dz \approx \sum_{q=0}^{n_q-1} \widetilde{w}_q f(\widetilde{\zeta}^{(q)})
$$

A particular quadrature rule $\mathcal{Q}=\{(\widetilde{\zeta}^{(q)},\widetilde{w}_q)\}_{q=0}^{n_q-1}$ is defined by the sets of points $\widetilde{\zeta}^{(i)}$ and corresponding weights $w_i$. Here we will consider Gauss-Legendre quadrature $\mathcal{Q}^{(\text{GL})}_{n_q}$, for which the points are the roots of the [Legendre polynomial](https://mathworld.wolfram.com/LegendrePolynomial.html) $P_{n_q}(\zeta)$ and the weights are given by $\widetilde{w}_q = \frac{2}{(1-(\widetilde{\zeta}^{(q)})^2)(P_{n_q}'(\widetilde{\zeta}^{(q)}))^2}$. The points and weights can be constructed with [numpy.polynomial.legendre.leggauss](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.legendre.leggauss.html):

```
points, weights = numpy.polynomial.legendre.leggauss(n)
```
The details of this construction are irrelevant for this course, but we need to have some understanding of how well the numerical scheme approximates the true value of the integral. Naturally one would expect that the quadrature approximates the integral better for larger numbers of points $n_q$. Crucially, Gauss-Legendre quadrature is exact if the function to be integrated is a polynomial of degree $2n_q-1$:

$$
\int_{-1}^{+1} p(z)\;dz = \sum_{q=0}^{n_q-1} \widetilde{w}_q p(\widetilde{\zeta}^{(q)})\qquad\text{for $p\in\mathcal{P}_{2n_q-1}$}
$$

We also call the degree of the highest polynomial that can be integrated exactly with a given quadrature rule the **degree of precision** or short "dop":
$$
\text{dop}(\mathcal{Q}^{(\text{GL})}_{n_q}) = 2n_q-1
$$
While so far we have only considered inetgration over the interval $[-1,+1]$, it turns out that integration over more general domains and higher-dimensional can be reduced to this case.

## Integration along a line
Next, imagine that we want to integrate a function along a straight line $\mathcal{C}\subset \mathbb{R}^2$ connecting two points $a,b\in \mathbb{R}^2$. To achieve this, pick a parametrisation $\gamma: [-1,1] \rightarrow \mathbb{R}^2$ of this line with $\gamma(-1)=a$, $\gamma(1)=b$

$$
\gamma(z) = \frac{1-z}{2}a+\frac{1+z}{2}b
$$

then

$$
\int_{\mathcal{C}} f(x)\;ds = \frac{\|b-a\|}{2} \int_{-1}^{+1} f(\gamma(z))\;dz
$$
Let $\mathcal{Q}_{n_q}^{(\text{GL})}=\{(\widetilde{\zeta}^{(q)},\widetilde{w}_q)\}_{q=0}^{n_q-1}$ be the Gauss-Legendre quadrature rule for the interval $[-1,+1]$. Then we obtain

$$
 \int_{\mathcal{C}} f(x)\;ds \approx \sum_{q=0}^{n_q-1} w_{\mathcal{C},q} f(\zeta_{\mathcal{C}}^{(q)})
$$

where the Gauss-Legendre quadrature rule on $\mathcal{C}$ is given by $\mathcal{Q}^{(\mathcal{C})}_{n_q}=\{(\zeta^{(q)}_{\mathcal{C}},w_{\mathcal{C},q})\}_{q=0}^{n_q-1}$ with

$$
\zeta_{\mathcal{C}}^{(q)} = \gamma(\widetilde{\zeta}^{(q)}) = \frac{1}{2}(1-\widetilde{\zeta}^{(q)})a + \frac{1}{2}(1+\widetilde{\zeta}^{(q)})b,\qquad
w_q = \|\gamma'(\zeta^{(q)})\| \widetilde{w}_{\mathcal{C},q} = \frac{\|b-a\|}{2} \widetilde{w}_q
$$
and
$$
\text{dop}(\mathcal{Q}^{(\mathcal{C})}_{n_q}) = \text{dop}(\mathcal{Q}^{(\text{GL})}_{n_q}) = 2n_q-1.
$$


## Two-dimensional quadrature for the reference triangle
To numerically integrate functions over the reference triangle $K$, first observe that $K$ is the image of the square $S=[-1,+1]\times [-1,+1]$ under the Duffy transform $\tau$ which maps a point $\widetilde{x}=(\widetilde{x}_0,\widetilde{x}_1)\in S$ to

$$
\begin{aligned}
\tau(\widetilde{x}) = \left(\frac{1}{2}(1+\widetilde{x}_0),\frac{1}{4}(1-\widetilde{x}_0)(1+\widetilde{x}_1)\right)^\top \in K
\end{aligned}
$$

### Integration over $\boldsymbol{S}$
Since $S=[-1,+1]\times[-1,+1]$ is the product of two intervals, we can integrate functions over $S$ by picking quadrature points and weights $\mathcal{Q}_{n_q}^{(S)}=\{(\widetilde{\zeta}^{(q)},\widetilde{w}_q)\}_{q=0}^{N_q-1}$ with $N_q = n_q(n_q+1)$ and

$$
\widetilde{\zeta}^{(q)} = \left(\widetilde{\zeta}^{(q_0)}_0,\widetilde{\zeta}^{(q_1)}_1\right)^\top,\quad \widetilde{w}_i = \widetilde{w}_{0,q_0}\cdot \widetilde{w}_{1,q_1} \qquad \text{where $q=n_q q_0+q_1$}.
$$

Here $\mathcal{Q}_{n_q+1}^{(\text{GL})} = \{(\widetilde{\zeta}^{(q_0)}_0,\widetilde{w}_{0,q_1})\}_{q_0=0}^{n_q}$ and $\mathcal{Q}_{n_q}^{(\text{GL})} =\{(\widetilde{\zeta}^{(q_1)}_0,\widetilde{w}_{1,q_1})\}_{q_1=0}^{n_q-1}$ are Gauss-Legendre quadrature rules with $n_q+1$ and $n_q$ points respectively (we need to integrate more accurately in the $0$-direction since an additional factor of $\widetilde{x}_0$ is introduced by the Duffy-transform). 

### Integration over $\boldsymbol{K}$
The quadrature rule $\mathcal{Q}_{n_q}^{(K)} = \{(\zeta^{(q)},w_q)\}_{q=0}^{N_q-1} = \tau(\mathcal{Q}^{(S)}_{n_q})$ over $K$ is then obtained as

$$
\begin{aligned}
\zeta^{(q)} &= \tau(\widetilde{\zeta}^{(q)}) = \left(\frac{1}{2}(1+\widetilde{\zeta}^{(q_0)}_0),\frac{1}{4}(1-\widetilde{\zeta}^{(q_0)}_0)(1+\widetilde{\zeta}^{(q_1)}_1)\right)^\top,\\
w_q &= \widetilde{w}_q \left|\det\left(\frac{\partial \tau}{\partial \widetilde{x}}\right)\right|_{\widetilde{x}=\widetilde{\zeta}^{(q)}} = \frac{1}{8}\widetilde{w}_{0,q_0}\widetilde{w}_{1,q_1}(1-\widetilde{\zeta}^{(q_0)}_0)
 \qquad \text{where $q=n_qq_0+q_1$.}
\end{aligned}
$$

The following figure shows the quadrature points on $S$ and $K$ for $n=2$ points.

![Quadrature points on $S$ and $K$](figures/quadrature.png)

## Implementation in Python

## FEM method on reference triangle
We can now implement a simple finite element method on the domain $\Omega=K$ defined by the reference triangle. For this, note that the entries of the stiffness matrix are given by:
$$
\begin{aligned}
A^{(h)}_{ij} = a(\phi_i,\phi_j) &= \int_K \left(\kappa \sum_{k=0}^{d-1}\frac{\partial\phi_i}{\partial x_k}(x) \frac{\partial\phi_j}{\partial x_k}(x) + \omega\; \phi_i(x) \phi_j(x)\right)\;dx\\
&\approx 
\sum_{q=0}^{N_q-1} w_q\left(\kappa \sum_{k=0}^{d-1}\frac{\partial\phi_i}{\partial x_k}(\zeta^{(q)}) \frac{\partial\phi_j}{\partial x_k}(\zeta^{(q)}) + \omega\; \phi_i(\zeta^{(q)}) \phi_j(\zeta^{(q)})\right)\\
&= \kappa \sum_{q=0}^{N_q-1}\sum_{k=0}^{d-1} w_q  T^\partial_{qik} (\boldsymbol{\zeta})T^\partial_{qjk} (\boldsymbol{\zeta})
+\omega \sum_{q=0}^{N_q-1} w_qT_{qi}(\boldsymbol{\zeta})T_{qj}(\boldsymbol{\zeta})
\end{aligned}
$$
Here $d=2$ is the dimension of the domain and $\{w_q,\zeta^{(q)}\}_{q=0}^{N_q-1}$ is a suitable quadrature rule on $K$. If we use a Lagrange finite element of polynomial degree $p$, then we need that make sure that the degree of precision of the quadrature rule is $2p$ to ensure that the product $\phi_i(x)\phi_j(x)$ is integrated exactly. Hence, we should use the quadrature rule $\mathcal{Q}_{p+1}^{(K)}$.

The entries of the right-hand side vector $\boldsymbol{b}^{(h)}$ are computed like this:
$$
\begin{aligned}
b^{(h)}_i = b(\phi_i) &= \int_K f(x)\phi_i(x)\;dx + \int_{\partial K} g(x)\phi_i(x)\;dx\\
&\approx \sum_{q=0}^{N_q-1} w_q f(\zeta^{(q)}) \phi_i(\zeta^{(q)}) + \sum_{\text{facets}\;F_m} \sum_{q=0}^{n_q-1 }w_{F_m,q} g(\zeta_{F_m}^{(q)})\phi_i(\zeta_{F_m}^{(q)}) \\
&= \sum_{q=0}^{N_q-1} w_q f_q(\boldsymbol{\zeta}) T_{qi}(\boldsymbol{\zeta}) + \sum_{\text{facets}\;F_m} \sum_{q=0}^{n_q-1 }w_{F_m,q} g_{q}(\boldsymbol{\zeta}_{F_m})T_{qi}(\boldsymbol{\zeta}_{F_m})
\end{aligned}
$$
with $f_q(\boldsymbol{\zeta}):=f(\zeta^{(q)})$ and $g_{q}(\boldsymbol{\zeta}_{F_m}) := g(\zeta_{F_m}^{(q)})$. We choose the quadrature rules $\mathcal{Q}_{n_q}^{(F_m)} = \{w_{F_m,q},\zeta^{(q)}_{F_m}\}_{q=0}^{n_q-1}$ with $n_q=p+1$ on the facets $F_m$.

The error $e^{(h)}(x)=u^{(h)}_{\text{exact}}(x)-u^{(h)}(x)$ is the difference between the exact and numerical solution. Furthermore, we have that

$$
e^{(h)}(x) = \sum_{j=0}^{d_p-1} \underbrace{\left(u^{(h)}_j - (u^{(h)}_{\text{exact}})_j\right)}_{=:e^{(h)}_j}\phi_j(x) = 
\sum_{j=0}^{d_p-1} e^{(h)}_j \phi_j(x)
$$
with $(u^{(h)}_{\text{exact}})_j = \ell_j(u_\text{exact})$ since the Lagrange element is nodal.

The square of the $L_2$ norm of the error is given by

$$
\begin{aligned}
\|e^{(h)}\|_{L_2}^2 &= \int_K (e^{(h)}(x))^2\;dx\\
&= \int_K\sum_{j,k=0}^{d_p-1} e^{(h)}_je^{(h)}_k \phi_j(x)\phi_k(x)\;dx\\
&\approx 
\sum_q ^{N_q-1}\sum_{j,k=0}^{d_p-1} w_q e^{(h)}_je^{(h)}_k \phi_j(\zeta^{(q)})\phi_k(\zeta^{(q)})\\
&= \sum_q ^{N_q-1}\sum_{j,k=0}^{d_p-1} w_q e^{(h)}_je^{(h)}_k T_{qj}T_{qk}.
\end{aligned}
$$

where $\mathcal{Q}_{n_q}^{(K)}=\{w_q\zeta^{(q)}\}_{q=0}^{N_q-1}$ is a suitable quadrature rule on $K$.

To test this, we use the method of manufactured solutions. For this, we pick a right-hand side $f(x)$ and boundary condition $g(x)$ such that the exact solution of $-\kappa \Delta u(x) + \omega u(x) = f(x)$ is given by

$$
u_{\text{exact}}(x) = \exp\left[-\frac{1}{2\sigma^2}(x-x_0)^2\right]
$$

A straightforward calculation shows that

$$
\begin{aligned}
f(x) &= \left(\left(2\frac{\kappa}{\sigma^2}+\omega\right)-\frac{\kappa}{\sigma^4}(x-x_0)^2\right) u_{\text{exact}}(x)
\\
g(x) &= -\frac{\kappa}{\sigma^2}n\cdot(x-x_0)u_{\text{exact}}(x)
\\
n &= \begin{cases}
\frac{1}{\sqrt{2}}\begin{pmatrix}1 \\ 1\end{pmatrix} & \text{for $x\in F_0$, i.e. $0\le x_0\le 1$, $x_0+x_1=1$}\\[3ex]
\begin{pmatrix}0 \\ -1\end{pmatrix} & \text{for $x\in F_1$, i.e. $0\le x_0\le 1$, $x_1=0$}\\[3ex]
\begin{pmatrix}-1 \\ 0\end{pmatrix} & \text{for $x\in F_2$, i.e. $x_0=0$, $0\le x_1\le 1$}
\end{cases}
\end{aligned}
$$

### Exercise
* Implement a method `assemble_lhs_reference_triangle(element, n_q)` which assembles the stiffness matrix $A^{(h)}$ using the Gauss-Legendre quadrature rule. The method should be passed:
  - An instance `element` of a subclass of `FiniteElement`
  - The number of points `n_q` used for the Gauss-Legendre quadrature
* Implement a method `assemble_rhs_reference_triangle(f, g, element, n_q)` which assembles the right-hand side vector $\boldsymbol{b}^{(h)}$ using the Gauss-Legenre quadrature rule. The method should be passed:
  - The function `f` which describes the right-hand side function $f(x)$
  - The function `g` which describes the Neumann boundary function $g(x)$
  - An instance `element` of a subclass of `FiniteElement`
  - The number of points `n_q` used for the Gauss-Legendre quadrature
* Implement a method `two_norm_reference_triangle(w, element, n_q)` which computes the $L_2$ norm of a function $w^{(h)}(x)=\sum_{j=0}^{d_p-1} w^{(h)}_j \phi_j(x)$. The method should be passed:
  - The vector $\boldsymbol{w}^{(h)}$ that defines the function $w^{(h)}(x)$ 
  - An instance `element` of a subclass of `FiniteElement`
  - The number of points `n_q` used for the Gauss-Legendre quadrature
* Solve the linear system $A^{(h)}\boldsymbol{u}^{(h)}=\boldsymbol{b}^{(h)}$ for the vector $\boldsymbol{u}^{(h)}$ by using [numpy.linalg.solve()](https://numpy.org/doc/2.0/reference/generated/numpy.linalg.solve.html)
* Apply the `tabulate_dofs()` method of the finite element class to the exact solution $u_{\text{exact}}(x)$ to obtain a vector $\boldsymbol{u}_{\text{exact}}^{(h)}$.
* Compute the error norm $\|e^{(h)}\|_{L_2}=\|\boldsymbol{u}_{\text{exact}}^{(h)}-\boldsymbol{u}^{(h)}\|_{L_2}$.
* How does $\|e^{(h)}\|_{L_2}$ depend on the polynomial degree $p$ of the Lagrange element?
* What happens for large values of $p$?