# Numerical quadrature
The weak form of the PDE is defined via suitable integrals such as $\int_\Omega v(x)f(x)\;dx$. In general, it is not possible to evaluate these integrals exactly. Furthermore, since the finite element discretisation (replacing $\mathcal{V}\mapsto \mathcal{V}_h$ and solving the associated linear algebra problem) already introduces an error, exact integration is not necessary, provided we can find an approximate integration method with errors that are of the same order of magnitude as the discretisation.

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

where the Gauss-Legendre quadrature rule on $\mathcal{C}$ is given by $\mathcal{Q}^{(\text{GL},\mathcal{C})}_{n_q}=\{(\zeta^{(q)}_{\mathcal{C}},w_{\mathcal{C},q})\}_{q=0}^{n_q-1}$ with

$$
\zeta_{\mathcal{C}}^{(q)} = \gamma(\widetilde{\zeta}^{(q)}) = \frac{1}{2}(1-\widetilde{\zeta}^{(q)})a + \frac{1}{2}(1+\widetilde{\zeta}^{(q)})b,\qquad
w_q = \|\gamma'(\zeta^{(q)})\| \widetilde{w}_{\mathcal{C},q} = \frac{\|b-a\|}{2} \widetilde{w}_q
$$
and
$$
\text{dop}(\mathcal{Q}^{(\text{GL},\mathcal{C})}_{n_q}) = \text{dop}(\mathcal{Q}^{(\text{GL})}_{n_q}) = 2n_q-1.
$$


## Two-dimensional quadrature for the reference triangle
To numerically integrate functions over the reference triangle $\widehat{K}$, first observe that $\widehat{K}$ is the image of the square $S=[-1,+1]\times [-1,+1]$ under the Duffy transform $\tau$ which maps a point $\widetilde{x}=(\widetilde{x}_0,\widetilde{x}_1)\in S$ to

$$
\begin{aligned}
\tau(\widetilde{x}) = \left(\frac{1}{2}(1+\widetilde{x}_0),\frac{1}{4}(1-\widetilde{x}_0)(1+\widetilde{x}_1)\right)^\top \in \widehat{K}
\end{aligned}
$$

(see figure below)

### Integration over $\boldsymbol{S}$
Since $S=[-1,+1]\times[-1,+1]$ is the product of two intervals, we can integrate functions over $S$ by picking quadrature points and weights $\mathcal{Q}_{n_q}^{(\text{GL},S)}=\{(\widetilde{\zeta}^{(q)},\widetilde{w}_q)\}_{q=0}^{N_q-1}$ with $N_q = n_q(n_q+1)$ and

$$
\widetilde{\zeta}^{(q)} = \left(\widetilde{\zeta}^{(q_0)}_0,\widetilde{\zeta}^{(q_1)}_1\right)^\top\in\mathbb{R}^2,\quad \widetilde{w}_i = \widetilde{w}_{0,q_0}\cdot \widetilde{w}_{1,q_1} \qquad \text{where $q=n_q q_0+q_1$}.
$$

Here $\mathcal{Q}_{n_q+1}^{(\text{GL})} = \{(\widetilde{\zeta}^{(q_0)}_0,\widetilde{w}_{0,q_1})\}_{q_0=0}^{n_q}$ and $\mathcal{Q}_{n_q}^{(\text{GL})} =\{(\widetilde{\zeta}^{(q_1)}_0,\widetilde{w}_{1,q_1})\}_{q_1=0}^{n_q-1}$ are Gauss-Legendre quadrature rules with $n_q+1$ and $n_q$ points respectively (we need to integrate more accurately in the $0$-direction since an additional factor of $\widetilde{x}_0$ is introduced by the Duffy-transform). 

### Integration over $\boldsymbol{\widehat{K}}$
The quadrature rule $\mathcal{Q}_{n_q}^{(\text{GL},\widehat{K})} = \{(\zeta^{(q)},w_q)\}_{q=0}^{N_q-1} = \tau(\mathcal{Q}^{(S)}_{n_q})$ over $\widehat{K}$ is then obtained as

$$
\begin{aligned}
\zeta^{(q)} &= \tau(\widetilde{\zeta}^{(q)}) = \left(\frac{1}{2}(1+\widetilde{\zeta}^{(q_0)}_0),\frac{1}{4}(1-\widetilde{\zeta}^{(q_0)}_0)(1+\widetilde{\zeta}^{(q_1)}_1)\right)^\top,\\
w_q &= \widetilde{w}_q \left|\det\left(\frac{\partial \tau}{\partial \widetilde{x}}\right)\right|_{\widetilde{x}=\widetilde{\zeta}^{(q)}} = \frac{1}{8}\widetilde{w}_{0,q_0}\widetilde{w}_{1,q_1}(1-\widetilde{\zeta}^{(q_0)}_0)
 \qquad \text{where $q=n_qq_0+q_1$.}
\end{aligned}
$$

The following figure shows the quadrature points on $S$ and $\widehat{K}$ for $n_q=2$.

![Quadrature points on $S$ and $\widehat{K}$](figures/quadrature.png)

Based on this construction we find that
$$
\text{dop}(\mathcal{Q}^{(\text{GL},\widehat{K})}_{n_q}) = \text{dop}(\mathcal{Q}^{(\text{GL})}_{n_q}) = 2n_q-1.
$$

## Implementation in Python

### Abstract base class
All quadrature rules are characterised by the weights and points. We therefore implement them as subclasses of an abstract base class `Quadrature` (in `fem/quadrature.py`) which has the following abstract properties:
* `nodes` the quadrature nodes $\{\xi^{(q)}\}_{q=0}^{n_q-1}$, represented by an array of shape $n_q\times 2$
* `weights` the quadrature weights $\{w_q\}_{q=0}^{n_q-1}$, represented by an array of length $n_q$
* `degree_of_precision` tegree of precision, i.e. the highest polynomial degree that can be integrated exactly

### Concrete implementations
The file `fem/quadrature.py` also contains specific subclasses

* A quadrature rule $\mathcal{Q}^{(\text{GL},\mathcal{C})}_{n_q}$ over line segments based on the Gauss-Legendre points can be implemented with `GaussLegendreQuadratureLineSegment(v_a, v_b, npoints)`. The following parameters are passed to the constructor:
    - `v_a` the start point $a$ of the line segment
    - `v_b` the end point $b$ of the line segment
    - `npoints` the number of points $n_q$
* A quadrature rule $\mathcal{Q}^{(\text{GL},\widehat{K})}_{n_q}$ over the reference triangle $\widehat{K}$ based on the Gauss-Legendre points can be implemented with `GaussLegendreQuadratureReferenceTriangle(npoints)`. The constructor is passed the number of points $n_q$.

## FEM method on reference triangle
We can now implement a simple finite element method on the domain $\Omega=\widehat{K}$ defined by the reference triangle.

### Stiffness matrix
For this, note that the entries of the stiffness matrix are given by:
$$
\begin{aligned}
A^{(h)}_{\ell k} = a(\phi_\ell,\phi_k) &= \int_{\widehat{K}} \left(\kappa \sum_{a=0}^{d-1}\frac{\partial\phi_\ell}{\partial x_a}(x) \frac{\partial\phi_k}{\partial x_a}(x) + \omega\; \phi_\ell(x) \phi_k(x)\right)\;dx\\
&\approx 
\sum_{q=0}^{N_q-1} w_q\left(\kappa \sum_{a=0}^{d-1}\frac{\partial\phi_\ell}{\partial x_a}(\zeta^{(q)}) \frac{\partial\phi_k}{\partial x_a}(\zeta^{(q)}) + \omega\; \phi_\ell(\zeta^{(q)}) \phi_k(\zeta^{(q)})\right)\\
&= \kappa \sum_{q=0}^{N_q-1}\sum_{a=0}^{d-1} w_q  T^\partial_{q\ell a} (\boldsymbol{\zeta})T^\partial_{qka} (\boldsymbol{\zeta})
+\omega \sum_{q=0}^{N_q-1} w_qT_{q\ell}(\boldsymbol{\zeta})T_{qk}(\boldsymbol{\zeta})
\end{aligned}
$$
Here $d=2$ is the dimension of the domain and $\{w_q,\zeta^{(q)}\}_{q=0}^{N_q-1}$ is a suitable quadrature rule on $\widehat{K}$. If we use a Lagrange finite element of polynomial degree $p$, then we need that make sure that the degree of precision of the quadrature rule is $2p$ to ensure that the product $\phi_i(x)\phi_j(x)$ is integrated exactly. Hence, we should use the quadrature rule $\mathcal{Q}_{p+1}^{(\text{GL},\widehat{K})}$.

### Right hand side vector
The entries of the right-hand side vector $\boldsymbol{b}^{(h)}$ are computed like this:
$$
\begin{aligned}
b^{(h)}_\ell = b(\phi_\ell) &= \int_{\widehat{K}} f(x)\phi_\ell(x)\;dx + \int_{\partial \widehat{K}} g(x)\phi_\ell(x)\;dx\\
&\approx \sum_{q=0}^{N_q-1} w_q f(\zeta^{(q)}) \phi_\ell(\zeta^{(q)}) + \sum_{\text{facets}\;F_m} \sum_{q=0}^{n_q-1 }w_{F_m,q} g(\zeta_{F_m}^{(q)})\phi_\ell(\zeta_{F_m}^{(q)}) \\
&= \sum_{q=0}^{N_q-1} w_q f_q(\boldsymbol{\zeta}) T_{q\ell}(\boldsymbol{\zeta}) + \sum_{\text{facets}\;F_m} \sum_{q=0}^{n_q-1 }w_{F_m,q} g_{q}(\boldsymbol{\zeta}_{F_m})T_{q\ell}(\boldsymbol{\zeta}_{F_m})
\end{aligned}
$$
with $f_q(\boldsymbol{\zeta}):=f(\zeta^{(q)})$ and $g_{q}(\boldsymbol{\zeta}_{F_m}) := g(\zeta_{F_m}^{(q)})$. We choose the quadrature rule $\mathcal{Q}_{n_q}^{(\text{GL},F_m)} = \{w_{F_m,q},\zeta^{(q)}_{F_m}\}_{q=0}^{n_q-1}$ with $n_q=p+1$ on the facets $F_m$.

### Error
The error $e^{(h)}(x)=u^{(h)}_{\text{exact}}(x)-u^{(h)}(x)$ is the difference between the exact and numerical solution. We can write $e^{(h)}$ as

$$
e^{(h)}(x) = u_{\text{exact}}(x) - \sum_{j=0}^{\nu-1} u^{(h)}_\ell \phi_\ell(x)
$$

The square of the $L_2$ norm of the error is given by

$$
\begin{aligned}
\|e^{(h)}\|_{L_2(\widehat{K})}^2 &= \int_{\widehat{K}} \left(u_{\text{exact}}(x) - \sum_{j=0}^{\nu-1} u^{(h)}_\ell \phi_\ell(x)\right)^2\;dx\\
&\approx 
\sum_q ^{N_q-1} w_q \left(u_{\text{exact}}(\zeta^{(q)}) - \sum_{\ell=0}^{\nu-1} u^{(h)}_\ell \phi_\ell(\zeta^{(q)})\right)^2 \\
&= \sum_q ^{N_q-1} w_q e_q^2\quad\text{with}\;\; e_q := u^{(\text{exact})}_q - \sum_{\ell=0}^{\nu-1} u^{(h)}_\ell T_{q\ell},\;u^{(\text{exact})}_q := u_{\text{exact}}(\zeta^{(q)}).
\end{aligned}\qquad(\dagger)
$$

where $\mathcal{Q}_{n_q}^{(\widehat{K})}=\{w_q,\zeta^{(q)}\}_{q=0}^{N_q-1}$ is a suitable quadrature rule on $\widehat{K}$.

### Numerical experiment
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
* Implement a method `assemble_lhs(element, n_q)` which assembles the stiffness matrix $A^{(h)}$ using the Gauss-Legendre quadrature rule. The method should be passed:
  - An instance `element` of a subclass of `FiniteElement`
  - The number of points `n_q` used for the Gauss-Legendre quadrature
* Implement a method `assemble_rhs(f, g, element, n_q)` which assembles the right-hand side vector $\boldsymbol{b}^{(h)}$ using the Gauss-Legenre quadrature rule. The method should be passed:
  - The function `f` which describes the right-hand side function $f(x)$
  - The function `g` which describes the Neumann boundary function $g(x)$
  - An instance `element` of a subclass of `FiniteElement`
  - The number of points `n_q` used for the Gauss-Legendre quadrature
* Implement a method `error_nrm(u, u_exact, element, n_q)` which computes the $L_2$ error norm $\|e^{(h)}\|_{L_2(\widehat{K})}$ by using the approximation in $(\dagger). The method should be passed:
  - The vector $\boldsymbol{u}^{(h)}$ that defines the function $u^{(h)}(x)$ 
  - A function $u_{\text{exact}}$ which represents the exact solution and which can be evaluated at arbirtrary points $\zeta\in \widehat{K}$
  - An instance `element` of a subclass of `FiniteElement`
  - The number of points `n_q` used for the Gauss-Legendre quadrature
* Solve the linear system $A^{(h)}\boldsymbol{u}^{(h)}=\boldsymbol{b}^{(h)}$ for the vector $\boldsymbol{u}^{(h)}$ by using [numpy.linalg.solve()](https://numpy.org/doc/2.0/reference/generated/numpy.linalg.solve.html)

#### Numerical experiments
* Apply the `tabulate_dofs()` method of the finite element class to the exact solution $u_{\text{exact}}(x)$ to obtain a vector $\boldsymbol{u}_{\text{exact}}^{(h)}$. For this pick the following parameters:
  - width of peak $\sigma = 0.5$
  - location of peak $x_0 = (0.6, 0.25)^\top$
  - Coefficient of diffusion term $\kappa = 0.9$
  - Coefficient of zero order term $\omega = 0.4$
* Compute the error norm $\|e^{(h)}\|_{L_2(\widehat{K})}$.
* How does $\|e^{(h)}\|_{L_2(\widehat{K})}$ depend on the polynomial degree $p$ of the Lagrange element?
* What happens for large values of $p$?

#### Practicalities

* Implement `assemble_lhs`, `assemble_rhs` and `error_nrm` in the file `algorithms.py`

You can use the Python functions given below. Note that the argument `x` can be a vector representing a single two-dimensional point or an array of shape $n\times 2$ which represents a collection of $n$ two-dimensional points.

```Python
def u_exact(x, sigma, x0):
    """Analytical solution

    :arg x: point at which the function is evaluated
    :arg sigma: width of peak
    :arg x0: location of peak"""
    return np.exp(
        -1 / (2 * sigma**2) * ((x[..., 0] - x0[0]) ** 2 + (x[..., 1] - x0[1]) ** 2)
    )


def f(x, kappa, omega, sigma, x0):
    """function to interpolate

    :arg x: point at which the function is evaluated
    :arg kappa: coefficient of diffusion term
    :arg omega: coefficient of zero-order term
    :arg sigma: width of peak
    :arg x0: location of peak
    """
    x_sq = (x[..., 0] - x0[0]) ** 2 + (x[..., 1] - x0[1]) ** 2
    return (2 * kappa / sigma**2 + omega - kappa / sigma**4 * x_sq) * u_exact(
        x, sigma, x0
    )


def g(x, kappa, sigma, x0):
    """boundary function

    :arg x: point at which the function is evaluated
    :arg kappa: coefficient of diffusion term
    :arg sigma: width of peak
    :arg x0: location of peak
    """
    if np.all(x[..., 1]) < 1e-12:
        # facet F_1
        n_dot_x = -(x[..., 1] - x0[1])
    elif np.all(x[..., 0]) < 1e-12:
        # facet F_2
        n_dot_x = -(x[..., 0] - x0[0])
    else:
        # facet F_0
        n_dot_x = (x[..., 0] - x0[0] + x[..., 1] - x0[1]) / np.sqrt(2)
    return -kappa / sigma**2 * n_dot_x * u_exact(x, sigma, x0)
```

To pass these functions to a method which expects a function $f(x)$ with a single argument $x$ you can use [`functools.partial`](https://docs.python.org/3/library/functools.html#functools.partial):

```Python
import functools
f_prime = functools.partial(x, kappa=0.9, omega=0.4, sigma=0.5, x0=[0.6, 0.25])
```

We can now call the function $f'(x)$ with a single argument $x$:
```Python
x = np.asarray([0.4,0.7])
f_prime(x)
```
