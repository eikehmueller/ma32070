<div align="center">
  <p style="font-size:32px;">MA32070 Exercises</p>
</div>

In the following we list some prior knowledge that is expected for this course together with some links to further information.

----

*&#169; Eike Mueller, University of Bath 2025. These notes are copyright of Eike Mueller, University of Bath. They are provided exclusively for educational purposes at the University and are to be downloaded or copied for your private study only. Further distribution, e.g. by upload to external repositories, is prohibited. html generated with [pandoc](https://pandoc.org/) using [easy-pandoc-templates](https://github.com/ryangrose/easy-pandoc-templates) under the [GPL-3.0.1 license](https://github.com/ryangrose/easy-pandoc-templates?tab=GPL-3.0-1-ov-file#readme)*

----

# Exercise: Cubic Lagrange element
Implement the cubic Lagrange element ($p=3$) in the class `CubicElement` by subclassing the abstract base class `FiniteElement`. The $\nu=10$ Lagrange points are in this order (see also figure above):

$$
\{\xi^{(j)}\}_{j=0}^{9}=
\left\{
\underbrace{
\begin{pmatrix}0\\[1ex]0\end{pmatrix},
\begin{pmatrix}1\\[1ex]0\end{pmatrix},
\begin{pmatrix}0\\[1ex]1\end{pmatrix}}_{\text{vertices}},
\underbrace{\begin{pmatrix}\frac{2}{3}\\[1ex]\frac{1}{3}\end{pmatrix},
\begin{pmatrix}\frac{1}{3}\\[1ex]\frac{2}{3}\end{pmatrix},
\begin{pmatrix}0\\[1ex]\frac{2}{3}\end{pmatrix},
\begin{pmatrix}0\\[1ex]\frac{1}{3}\end{pmatrix},
\begin{pmatrix}\frac{1}{3}\\[1ex]0\end{pmatrix},
\begin{pmatrix}\frac{2}{3}\\[1ex]0\end{pmatrix}}_{\text{facets}},
\underbrace{\begin{pmatrix}\frac{1}{3}\\[1ex]\frac{1}{3}\end{pmatrix}}_{\text{interior}}\right\}
$$

You can use the following 10 monomials:

$$
\{\theta_\ell(x)\}_{\ell=0}^{9} = \{1,x_0,x_1,x_0^2,x_0x_1,x_1^2,x_0^3,x_0^2x_1,x_0x_1^2,x_1^3\}
$$

* Your class should store the Lagrange points in an attribute `_nodal_points`
* Your class should contain a method `vandermonde_matrix(zeta,grad=False)` which accepts as an argument a $n\times 2$ matrix of $n$ two-dimensional points. The method should compute the $n\times \nu$ matrix $V(\boldsymbol{\zeta})$ if `grad=False` and the $n\times \nu\times 2$ tensor $V^\partial(\boldsymbol{\zeta})$ if `grad=True`. If only a single point is passed to the method, it should return a vector of length $\nu$ for `grad=False` and a $\nu\times 2$ matrix for `grad=True`.
* Use the `vandermonde_matrix()` method together with `_nodal_points` to construct the coefficient matrix `C`
* Use the coefficient matrix `C` and the `vandermonde_matrix()` method to tabulate the basis functions and their gradients by using the expressions above. You might find the [`numpy.einsum()`](https://numpy.org/doc/2.2/reference/generated/numpy.einsum.html) method useful to compute $T^\partial(\boldsymbol{\zeta})$

## Practicalities
* Save your implementation in the file `cubicelement.py`
* Use [pytest](https://docs.pytest.org/) to develop a suite of suitable tests to check that your implementation is correct. For this, create a file `test_cubicelement.py` in the same directory

# Exercise: Local assembly on the reference triangle
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

## Numerical experiments
* Apply the `tabulate_dofs()` method of the finite element class to the exact solution $u_{\text{exact}}(x)$ to obtain a vector $\boldsymbol{u}_{\text{exact}}^{(h)}$. For this pick the following parameters:
  - width of peak $\sigma = 0.5$
  - location of peak $x_0 = (0.6, 0.25)^\top$
  - Coefficient of diffusion term $\kappa = 0.9$
  - Coefficient of zero order term $\omega = 0.4$
* Compute the error norm $\|e^{(h)}\|_{L_2(\widehat{K})}$.
* How does $\|e^{(h)}\|_{L_2(\widehat{K})}$ depend on the polynomial degree $p$ of the Lagrange element?
* What happens for large values of $p$?

## Practicalities

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

# Exercise: PETSc sparse matrices
Create two $3\times 3$ sparse PETSc matrices $A$, $B$.

By using suitable functions (see [`petsc4py.Mat` documentation](https://petsc.org/release/petsc4py/reference/petsc4py.PETSc.Mat.html)), compute
* $AB$ 
* $AB^\top$
* $A+B$

# Exercise: PETSc solver options
Instead of $P=D$, we could also use the lower triangular part of $A$ and set $P=D+L$ where
$$
L = \begin{cases}
A_{ij} & \text{if $i<j$} \\
0 & \text{otherwise}
\end{cases}
$$
Convince yourself that for a given vector $\boldsymbol{r}$ the equation $(D+L)z=r$ can be solved row-by-row, i.e. by computing first $\boldsymbol{z}_0 = \boldsymbol{r}_0/A_{00}$, then computing $\boldsymbol{z}_1 = (\boldsymbol{r}_1 - A_{10}\boldsymbol{z}_0)/A_{11}$, $\boldsymbol{z}_2=(\boldsymbol{r}_2 - A_{20}\boldsymbol{z}_0 - A_{21}\boldsymbol{z}_1)/A_{22}$ and so on. The corresponding preconditioner is also known as the successive overrelaxation (SOR) method. It can be chosen by setting `-pc_type sor``. Run the code with this preconditioner - how does the number of iterations change?