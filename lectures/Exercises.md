<div align="center">
  <p style="font-size:32px;">MA32070 Exercises</p>
</div>

----

*&#169; Eike Mueller, University of Bath 2025. These notes are copyright of Eike Mueller, University of Bath. They are provided exclusively for educational purposes at the University and are to be downloaded or copied for your private study only. Further distribution, e.g. by upload to external repositories, is prohibited. html generated with [pandoc](https://pandoc.org/) using [easy-pandoc-templates](https://github.com/ryangrose/easy-pandoc-templates) under the [GPL-3.0.1 license](https://github.com/ryangrose/easy-pandoc-templates?tab=GPL-3.0-1-ov-file#readme)*

----

# Exercise 1: Cubic Lagrange element 

#### Set: week 2
#### Due: at the end of week 4

## Task
Implement the cubic Lagrange element ($p=3$) in the class `CubicElement` by subclassing the abstract base class `FiniteElement`. The $\nu=10$ Lagrange points are in this order (see also figure below):

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

![Lagrange points for cubic element](figures/lagrange_nodes_cubic.svg)

The corresponding nodal basis functions with $\phi_\ell(\xi^{(k)}) = \delta_{jk}$ are
$$
\begin{aligned}
\phi_0(x) &= 1 - \frac{11}{2} (x_0 + x_1) + 9 (x_0 + x_1)^2 - \frac{9}{2} (x_0 + x_1)^3,\\
\phi_1(x) &= x_0\left(1 - \frac{9}{2} x_0 + \frac{9}{2} x_0^2 \right),\\
\phi_2(x) &= x_1\left(1 - \frac{9}{2} x_1 + \frac{9}{2} x_1^2\right),\\
\phi_3(x) &= -\frac{9}{2} x_0 x_1 (1 - 3 x_0),\\
\phi_4(x) &= -\frac{9}{2} x_0 x_1  (1 - 3 x_1),\\
\phi_5(x) &= -\frac{9}{2} x_1 (1 - x_0 - 4 x_1 + 3 x_1 * (x_0 + x_1)),\\
\phi_6(x) &= \frac{9}{2} x_1 (2 - 5 (x_0 + x_1) + 3 (x_0 + x_1)^2),\\
\phi_7(x) &= \frac{9}{2} x_0 (2 - 5 (x_0 + x_1) + 3 (x_0 + x_1)^2),\\
\phi_8(x) &= -\frac{9}{2} x_0 (1 - x_1 - 4 x_0 + 3 x_0 (x_0 + x_1)),\\
\phi_9(x) &= 27 x_0 x_1 (1 - x_0 - x_1)
\end{aligned}
$$

You can use the following 10 monomials for the construction of the Vandermonde matrix:

$$
\{\theta_\ell(x)\}_{\ell=0}^{9} = \{1,x_0,x_1,x_0^2,x_0x_1,x_1^2,x_0^3,x_0^2x_1,x_0x_1^2,x_1^3\}
$$
## Requirements
* Your class should store the Lagrange points in an attribute `_nodal_points`, this can be done for example in the `__init__()` method.
* Your class should contain a method `vandermonde_matrix(zeta,grad=False)` which accepts as an argument a $n\times 2$ array of $n$ two-dimensional points. The method should compute the $n\times \nu$ matrix $V(\boldsymbol{\zeta})$ defined in the lectures if `grad=False` and the $n\times \nu\times 2$ tensor $V^\partial(\boldsymbol{\zeta})$ if `grad=True`. If only a single point is passed to the method, it should return a vector of length $\nu$ for `grad=False` and an $\nu\times 2$ array for `grad=True`.
* In the `__init__()` method, use the `vandermonde_matrix()` method together with `_nodal_points` to construct the coefficient matrix $C$ defined in the lectures. Store this matrix in a class attribute `_coefficients`
* Use `_coefficients` and the `vandermonde_matrix()` method to tabulate the basis functions and their gradients to implement the `tabulate()` and `tabulate_gradient()` methods. You might find the [`numpy.einsum()`](https://numpy.org/doc/2.2/reference/generated/numpy.einsum.html) method useful to compute $T^\partial(\boldsymbol{\zeta})$
* Verify that the methods `dofmap()` and `inverse_dofmap()`, which are inherited from the abstract base class behave as expected. For this, inspect the output of `inverse_dofmap(ell)` for $\ell=0,1,2,\dots,9$
* Make sure you implement all other methods/properties that are marked as abstract in the base class.
  
Your code should pass the tests below, which verify correctness for special cases. Use [pytest](https://docs.pytest.org/) to add further tests to verify that your implementation is correct. In particular, you should check that
* `tabulate()` correctly computes $\phi_\ell(\xi^{(k)}) = \delta_{\ell k}$ where $\xi^{(k)}$ are the nodal points.
* `tabulate_dofs()` correctly computes $\lambda_\ell(\phi_k) = \delta_{\ell k}$
  
## Practicalities
* Save your implementation in the file `cubicelement.py` and the tests in `test_cubicelement.py` in the same directory `ma32070/exercise1`
* Zip this directory which contains `cubicelement.py` and `test_cubicelement.py`. For this, change to `ma32070/` and run `tar czvf exercise1.tgz exercise1`
* Upload the resulting file `exercise1.tgz` to the submission point on moodle

## Tests
```python
def test_ndof_per_vertex():
    """Check that the number of unknowns per vertex is set correctly"""
    element = CubicElement()
    assert element.ndof_per_vertex == 1


def test_ndof_per_facet():
    """Check that the number of unknowns per facet is set correctly"""
    element = CubicElement()
    assert element.ndof_per_facet == 2


def test_ndof_per_interior():
    """Check that the number of unknowns in the interior is set correctly"""
    element = CubicElement()
    assert element.ndof_per_interior == 1


def test_tabulate_dofs():
    """Check that dof-tabulation of exp(-x)*(2+sin(x)) is correct"""
    element = CubicElement()
    expected = [
        2.00000000, 0.73575888, 2.84147098, 1.19482160, 1.87614395,
        2.61836980, 2.32719470, 1.43306262, 1.02683424, 1.66750787,
    ]
    tabulated = element.tabulate_dofs(
        lambda x: np.exp(-x[..., 0]) * (2 + np.sin(x[..., 1]))
    )
    assert np.allclose(expected, tabulated, rtol=1.0e-8)


def test_tabulate_single_point():
    """Check that tabulation of all dofs at a single point is correct"""
    element = CubicElement()
    zeta = [0.18, 0.43]
    expected = [
        -0.0275145, 0.060444, -0.0442685, -0.160218, 0.101007,
        0.2188485, 0.1282905, 0.053703, -0.145314, 0.815022,
    ]
    tabulated = element.tabulate(zeta)
    assert np.allclose(tabulated, expected, rtol=1.0e-12)


def test_tabulate_multiple_points():
    """Check that tabulation of all dofs at multiple points is correct"""
    element = CubicElement()
    zeta = [[0.18, 0.43], [0.72, 0.21], [0.4, 0.31]]
    tabulated = element.tabulate(zeta)
    expected = [
        [
            -0.0275145, 0.060444, -0.0442685, -0.160218, 0.101007,
            0.2188485, 0.1282905, 0.053703, -0.145314, 0.815022,
        ],
        [
            0.0494935, 0.066816, 0.0532245, 0.789264, -0.251748,
            -0.0244755, -0.0522585, -0.179172, 0.263088, 0.285768,
        ],
        [
            0.0213005, -0.032, 0.0116095, 0.1116, -0.03906,
            -0.0283185, -0.0525915, -0.06786, 0.1044, 0.97092,
        ],
    ]
    assert np.allclose(tabulated, expected, rtol=1.0e-12)


def test_tabulate_gradient_single_point():
    """Check that tabulation of gradient of all dofs at a single point is correct"""
    element = CubicElement()
    zeta = [0.18, 0.43]
    tabulated = element.tabulate_gradient(zeta)
    expected = [
        [0.45665, 0.45665], [-0.1826, 0.0],
        [0.0, -0.37385], [0.1548, -0.3726],
        [0.56115, 1.2798], [-0.56115, 2.21175],
        [-2.5929, -2.29455], [-0.78705, -1.0854],
        [0.513, 0.3726], [2.4381, -0.1944],
    ]
    assert np.allclose(tabulated, expected, rtol=1.0e-12)


def test_tabulate_gradient_multiple_points():
    """Check that tabulation of gradient of all dofs at multiple points is correct"""
    element = CubicElement()
    zeta = [[0.18, 0.43], [0.72, 0.21], [0.4, 0.31]]
    tabulated = element.tabulate_gradient(zeta)
    expected = [
        [
            [0.45665, 0.45665], [-0.1826, 0.0],
            [0.0, -0.37385], [0.1548, -0.3726],
            [0.56115, 1.2798], [-0.56115, 2.21175],
            [-2.5929, -2.29455], [-0.78705, -1.0854],
            [0.513, 0.3726], [2.4381, -0.1944],
        ],
        [
            [-0.43615, -0.43615], [1.5184, 0.0],
            [-0.0, -0.29465], [3.1374, 3.7584],
            [-0.34965, 0.8424], [0.34965, 0.43155],
            [0.5481, 0.29925], [1.63035, 1.8792],
            [-2.7126, -3.7584], [-3.6855, -2.7216],
        ],
        [
            [0.47465, 0.47465], [-0.44, 0.0],
            [0.0, -0.49265], [1.953, 0.36],
            [-0.09765, 1.548], [0.09765, 1.21995],
            [-1.0323, -1.20195], [-1.50165, -1.332],
            [1.467, -0.36], [-0.9207, -0.216],
        ],
    ]
    assert np.allclose(tabulated, expected, rtol=1.0e-12)
```

# Exercise 2: Local assembly on the reference triangle

#### Set: week 3
#### Due: end of week 5

## Tasks
### Implementation
* Implement a method `assemble_lhs(element, n_q)` which assembles the stiffness matrix $A^{(h)}$ using the Gauss-Legendre quadrature rule. The method should be passed:
  - An instance `element` of a subclass of `FiniteElement`
  - The number of points `n_q` used for the Gauss-Legendre quadrature
* Implement a method `assemble_rhs(f, g, element, n_q)` which assembles the right-hand side vector $\boldsymbol{b}^{(h)}$ using the Gauss-Legendre quadrature rule. The method should be passed:
  - The function `f` which describes the right-hand side function $f(x)$
  - The function `g` which describes the Neumann boundary function $g(x)$
  - An instance `element` of a subclass of `FiniteElement`
  - The number of points `n_q` used for the Gauss-Legendre quadrature
* Implement a method `error_nrm(u_numerical, u_exact, element, n_q)` which computes the $L_2$ error norm $\|e_h\|_{L_2(\widehat{K})}$ by using the approximation from the lecture. The method should be passed:
  - The vector `u_numerical` which is the dof-vector $\boldsymbol{u}^{(h)}$ that represents the function $u_h(x)$ 
  - A function `u_exact` which represents the exact solution $u_{\text{exact}}(x)$ and which can be evaluated at arbitrary points $\zeta\in \widehat{K}$
  - An instance `element` of a subclass of `FiniteElement`
  - The number of points `n_q` used for the Gauss-Legendre quadrature

Use the methods `assemble_lhs()`, `assemble_lhs()` and `error_nrm()` in a main program which

* assembles the matrix $A^{(h)}$ and right hand $\boldsymbol{b}^{(h)}$,
* solves the linear system $A^{(h)}\boldsymbol{u}^{(h)}=\boldsymbol{b}^{(h)}$ for the vector $\boldsymbol{u}^{(h)}$ by using [`numpy.linalg.solve()`](https://numpy.org/doc/2.0/reference/generated/numpy.linalg.solve.html)
* computes the error norm $\|e_h\|_{L_2(\widehat{K})}$,
* visualises the solution and error.

### Numerical experiments

Study the accuracy of the solution for different polynomial degrees $p$. As the exact solution, pick a Gaussian
$$
u_{\text{exact}}(x) = \exp\left[-\frac{1}{2\sigma^2}\left\|\boldsymbol{x}-\boldsymbol{x}_0\right\|_2^2\right]
$$
with a peak of width $\sigma = 0.5$ centred at $x_0 = (0.6, 0.25)^\top$.

In the weak form which defines the PDE $-\kappa \Delta u + \omega u = f$, use the parameters $\kappa = 0.9$ and $\omega = 0.4$.

Compute the error norm $\|e_h\|_{L_2(\widehat{K})}$ and visualise the solution and error for polynomial degrees $p=1$ and $p=3$.

## Practicalities

* Implement the methods `assemble_lhs()`, `assemble_rhs()` and `error_nrm()` in the file `algorithms.py` in the directory `ma32070/exercise2`
* Use these methods in the file `driver.py` in the same directory
* Create a file `solution.pdf` file with the following content in the same directory:
    - a brief description of how you implemented and tested your code
    - a table which lists $\|e_h\|_{L_2(\widehat{K})}$ for $p=1$ and $p=3$
    - plots of the solution and error for $p=1$ and $p=3$
* Zip the directory which contains `algorithms.py`, `driver.py` and `solution.pdf`. For this, change to `ma32070/` and run `tar czvf exercise2.tgz exercise2`
* Upload the resulting file `exercise2.tgz` to the submission point on moodle

## Hints
* You can import the `LinearElement` and quadrature rules provided in the finite element library with 
 
```Python
from fem.linearelement import LinearElement
from fem.quadrature import (
    GaussLegendreQuadratureLineSegment,
    GaussLegendreQuadratureReferenceTriangle,
)
```
* The method `plot_solution(u_numerical, u_exact, element, filename)` in `fem.utilities` can be used to visualise the solution and the error; the result is written to a file. Look at the documentation of the method to understand how it is used.
* You can use the Python functions $u_{\text{exact}}(x)$, $f(x)$ and $g(x)$ given below. Note that the argument `x` can be a vector representing a single two-dimensional point or an array of shape $n\times 2$ which represents a collection of $n$ two-dimensional points.

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
    """Function f(x) for right hand side

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
    """Boundary function g(x)

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
f_star = functools.partial(f, kappa=0.9, omega=0.4, sigma=0.5, x0=[0.6, 0.25])
```

We can now call the function $f^*(x)$ with a single argument $x$:
```Python
x = np.asarray([0.4,0.7])
f_star(x)
```

# Bonus Exercise: Three point quadrature (not marked)

#### Set: week 3

## Task
Consider the following three-point quadrature on the reference triangle $\widehat{K}$:

$$
\begin{aligned}
w_0 &= \frac{1}{6}, & w_1 &= \frac{1}{6}, & w_2 &= \frac{1}{6}\\
\xi^{(0)} &= \begin{pmatrix}\frac{1}{6} \\[1.5ex] \frac{1}{6}\end{pmatrix}, &
\xi^{(1)} &= \begin{pmatrix}\frac{2}{3} \\[1.5ex] \frac{1}{6}\end{pmatrix}, &
\xi^{(2)} &= \begin{pmatrix}\frac{1}{6} \\[1.5ex] \frac{2}{3}\end{pmatrix}
\end{aligned}
$$

The degree of precision of this rule is $2$, i.e. it is exact for polynomials of the form
$$
p(x_0,x_1) = \sum_{\substack{s_0,s_1\\s_0+s_1\le 2}} a_{s_0,s_1} x_0^{s_0} x_1^{s_1}
$$

### Implementation
Implement this quadrature rule in the class `ThreePointQuadratureReferenceTriangle`, which should be a subclass of the abstract base class `Quadrature`. Write a suitable test which verifies that the implemention is correct. For this, observe that

$$
\int_{\widehat{K}} x_0^{s_0} x_1^{s_1}\;dx_0\;dx_1 = \frac{s_0!s_1!}{(s_0+s_1+2)!}
$$

Use the [`@pytest.mark.parametrize`](https://docs.pytest.org/en/stable/example/parametrize.html) decorator to write a suitable, parametrised test which verifies this:

```python
import pytest

@pytest.mark.parametrize(
    "s, expected",
    [
        [[0, 0], 1 / 2],
        # Add more cases here
    ],
)
def test_threepoint_quadrature_monomial(s, expected):
    """Check that three point quadrature is exact for monomial x_0^{s_0} x_1^{s_1}
    
    :arg s: (s_0,s_1) = powers of x_0,x_1
    :arg expected: exact result s_0! s_1! / (s_0+s_1+2)!
    """
    # Add your own test code here
```

### Theory
Confirm that the quadrature rule is exact for polynomials of degree less than or equal 2, i.e. that

$$
\sum_{q=0}^{2} w_q (\xi_0^{(q)})^{s_0}(\xi_1^{(q)})^{s_1}= \int_{\widehat{K}} x_0^{s_0} x_1^{s_1}\;dx_0\;dx_1\qquad \text{for $s_0,s_1\ge 0$ and $s_0+s_1\le 2$}.
$$

## Practicalities
* Implement `ThreePointQuadratureReferenceTriangle` in a file `threepointquadrature.py` in the directory `ma32070/bonus_exercise`
* Implement the tests in `test_threepointquadrature.py` in the same directory

# Exercise 3: Computation of global $L_2$-error

#### Set: week 6
#### Due: end of week 7

## Task
As for the simplified case where $\Omega=\widehat{K}$ is the reference triangle, the error $e_h(x)=u_{\text{exact}}(x)-u_h(x)$ is the difference between the exact solution and numerical solution $u_h(x)$. Expanding $u_h(x)$ in terms of the basis functions $\Phi_{\ell_{\text{globa;}}}(x)$, we can write the error $e_h$ as

$$
e_h(x) = u_{\text{exact}}(x) - \sum_{\ell_{\text{global}}=0}^{n-1} u^{(h)}_{\ell_{\text{global}}} \Phi^{(h)}_{\ell_{\text{global}}}(x).
$$

The square of the $L_2$ norm of the error can be computed by summing over all triangles in the mesh

$$
\begin{aligned}
\|e_h\|_{L_2(\Omega)}^2 &= \int_{\Omega} \left(u_{\text{exact}}(x) - \sum_{\ell_{\text{global}}=0}^{n-1} u^{(h)}_{\ell_{\text{global}}} \phi_{\ell_{\text{global}}}(x)\right)^2\;dx\\
&= \sum_{K\in\Omega_h} \int_{K} \left(u_{\text{exact}}(x) - \sum_{\ell=0}^{\nu-1} u^{(h)}_{\ell_{\text{global}}} \Phi^{(h)}_{\ell_{\text{global}}}(x)\right)^2\;dx\\
\end{aligned}
$$

Changing variables to integrate over the reference cell $\widehat{K}$ this leads to

$$
\begin{aligned}
\|e_h\|_{L_2(\Omega)}^2 &= \sum_{K\in\Omega_h} \int_{\widehat{K}} \left(\widehat{u}_{K,\text{exact}}(\widehat{x}) - \sum_{\ell=0}^{\nu-1} u^{(h)}_{\ell_{\text{global}}} \phi_{\ell}(\widehat{x})\right)^2\left|\det{J(\widehat{x})}\right|\;dx
\end{aligned}
$$

with $\ell_{\text{global}}=\ell_{\text{global}}(\alpha,\ell)$ the global index corresponding to the local dof-index $\ell$ in the cell with index $\alpha$ and $\widehat{u}_{K,\text{exact}} = u_{\text{exact}}\circ X_K$ the pullback of the exact solution to the cell $K$. 

Finally, we approximate integration by numerical quadrature to obtain

$$
\begin{aligned}
\|e_h\|_{L_2(\Omega)}^2 &\approx 
\sum_{K\in\Omega_h}\sum_{q=0} ^{N_q-1} w_q \left(\widehat{u}_{K,\text{exact}}(\zeta^{(q)}) - \sum_{\ell=0}^{\nu-1} u^{(h)}_{\ell_{\text{global}}} \phi_\ell(\zeta^{(q)})\right)^2 \left|\det{J(\zeta^{(q)})}\right|.
\end{aligned}
$$

where $\mathcal{Q}_{n_q}^{(\widehat{K})}=\{w_q,\zeta^{(q)}\}_{q=0}^{N_q-1}$ is a suitable quadrature rule on $\widehat{K}$.

This leads to the following procedure:

### Algorithm: Computation of global $L_2$ error
1. Initialise $S \gets 0$
1. For all cells $K$ **do**:
1. $~~~~$ Extract the coordinate dof-vector $\overline{\boldsymbol{X}}$ with $\overline{X}_{\ell^\times} = X_{\ell^\times_\text{global}(\alpha,{\ell^\times})}$ where $\alpha$ is the index of cell $K$
1. $~~~~$ Extract the local dof-vector $\overline{\boldsymbol{u}}$ with $\overline{u}_{\ell} = u^{(h)}_{\ell_\text{global}(\alpha,\ell)}$
1. $~~~~$ For all quadrature points $q$ **do**:
1. $~~~~~~~~$ Compute the determinant $D_q$ of the Jacobian $J(\xi^{(q)})$ with $J_{ab}(\xi^{(q)}) = \sum_{\ell^\times} \overline{X}_{\ell^\times} T^{\times\partial}_{q\ell^\times ab}$
1. $~~~~~~~~$ Compute $(x_K^{(q)})_a = \sum_{\ell^\times} T^\times_{q\ell^\times a} \overline{X}_{\ell^\times}$ and evaluate $u^{\text{(exact)}}_q = \widehat{u}_{K,\text{exact}}(\xi^{(q)}) = u_{\text{exact}}(x_K^{(q)})$
1. $~~~~~~~~$ Compute $e_q = u^{\text{(exact)}}_q - \sum_{\ell=0}^{\nu-1}T_{q\ell} \overline{u}_\ell$
2. $~~~~~~~~$ Update $S \gets S + w_q e_q^2 D_q$
3. $~~~~$ **end do**
4.  **end do**

### Implementation
Write a method `error_nrm(u_h, u_exact, quad)` which implements the above algorithm. Your method should take the following parameters:
* An object `u_h` of class `Function` which represents the numerical solution $u^{(h)}$
* A Python function `u_hexact` which can represents the exact solution $u_{\text{exact}}(x)$ and which be evaluated at arbitrary points $x\in \Omega$ 
* A suitable quadrature rule `quad`

You can use the following main program, which solves the PDE $-\kappa\Delta u + \omega u=f$ for $\kappa=0.9$, $\omega=0.4$ with the right hand side $f(x)$ chosen such that the exact solution is

$$
u_{\text{exact}}(x) = \cos(2\pi x_0)\cos(4\pi x_1)
$$

```Python
"""Solve finite element model problem and compute global L2 error"""

import numpy as np

from fem.utilitymeshes import RectangleMesh
from fem.linearelement import LinearElement
from fem.utilities import measure_time
from fem.functionspace import FunctionSpace
from fem.function import Function, CoFunction
from fem.algorithms import assemble_rhs, assemble_lhs
from fem.quadrature import GaussLegendreQuadratureReferenceTriangle
from algorithms import error_nrm

def f(x):
    """Right hand side

    :arg x: point at which to evaluate the function
    """
    return (
        ((2**2 + 4**2) * np.pi**2 * kappa + omega)
        * np.cos(2 * np.pi * x[0])
        * np.cos(4 * np.pi * x[1])
    )


def u_exact(x):
    """Exact solution

    :arg x: point at which to evaluate the function
    """
    return np.cos(2 * np.pi * x[0]) * np.cos(4 * np.pi * x[1])


# Number of mesh refinements
nref = 5
# Coeffcient of diffusion term
kappa = 0.9
# Coefficient of zero order term
omega = 0.4

# Finite element
element = LinearElement()

# Mesh
mesh = RectangleMesh(Lx=1, Ly=1, nref=nref)
# Function space
fs = FunctionSpace(mesh, element)
print(f"nref = {nref}")
print(f"grid spacing h = {np.sqrt(2)/2**nref}")
print(f"number of unknowns = {fs.ndof}")
# Quadrature rule
quad = GaussLegendreQuadratureReferenceTriangle(2)

# Construct right hand side
# (s_0^2 + s_1^2)*pi^2*u(x)
b_h = CoFunction(fs)
with measure_time("assemble_rhs()"):
    assemble_rhs(f, b_h, quad)

# Numerical solution
u_h = Function(fs, "u_numerical")

# Stiffness matrix
with measure_time("assemble_lhs()"):
    stiffness_matrix = assemble_lhs(fs, quad, kappa, omega)

# Solve linear system A^{(h)} u^{(h)} = b^{(h)}
with measure_time("solve()"):
    u_h.data[:] = np.linalg.solve(stiffness_matrix, b_h.data)

error_norm = error_nrm(u_h, u_exact, quad)

print()
print(f"error = {error_norm}")
```

### Numerical experiments
Run your code for different problem sizes by setting $n_{\text{ref}} = 3,4,5,6,7$.

#### Runtime measurements
Produce a table which shows the time spent in the following parts of the code
* Assembly of stiffness matrix $A^{(h)}$ in `assemble_lhs()`
* Assembly of right-hand side $\boldsymbol{b}^{(h)}$ in `assemble_rhs()`
* Solution of the linear system $A^{(h)}\boldsymbol{u}^{(h)}=\boldsymbol{b}^{(h)}$
  
as a function of the number of unknowns $n_{\text{dof}}$, which can be obtained as `fs.ndof` if `fs` is an instance of `FunctionSpace`. In the above main program, the `measure_time` decorator is used to time specific components of the code.

How does the time spent in different parts of the code increase as a function of $n_{\text{dof}}$?

#### Convergence
Visualise the norm $\|e_h\|_{L_2(\Omega)}$ of the $L_2$ error as a function of the grid spacing $h=\sqrt{2}\cdot 2^{-n_{\text{ref}}}$ in a log-log plot. 

Repeat this experiment for the `CubicElement` and $n_{\text{ref}}=3,4,5,6$; remember to adapt the order of the quadrature appropriately.

Assuming that for $h\ll 1$ we have
$$
\|e_h\|_{L_2(\Omega)}\approx C h^{\alpha}
$$
which empirical rate of convergence $\alpha$ do you observe in the two cases?

## Practicalities
* Save your implementation of `error_nrm()`in the file `algorithms.py` and the main program (copied from above) in `driver.py` in the same directory `ma32070/exercise3`
* Write a brief report (no more than 2 pages) which describes your implementation and presents and discusses your numerical results. Include this as a single `.pdf` file called `solution.pdf` in the same directory `ma32070/exercise3`
* Zip this directory which contains `algorithms.py`, `driver.py` and `solution.pdf`. For this, change to `ma32070/` and run `tar czvf exercise3.tgz exercise3`
* Upload the resulting file `exercise3.tgz` to the submission point on moodle
   
# Exercise 4: Computational cost of backsubstitution

#### Set: week 7
#### Due: end of week 8

## Solution of upper triangular systems
Consider a $n\times n$ upper triangular system, such as in the following example ($n=5$):
$$
\begin{aligned}
A &=  \begin{pmatrix}
1.096 & 0.3391 & 0.0632 & 0.0555 & 0.2176\\
\textcolor{lightgray}{0} & 1.197 & 0.06495 & 0.3045 & 0.01172\\
\textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 1.127 & 0.0008768 & 0.1785\\
\textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 1.263 & 0.3005\\
\textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 1.047
\end{pmatrix} &
    b &=  \begin{pmatrix}
0.23\\
0.4326\\
0.1808\\
0.7749\\
0.7747
\end{pmatrix}
\end{aligned}
$$

More generally, an upper triangular matrix satisfies $A_{ij} = 0$ for all $i>j$; we also assume that $A_{ii}\neq 0$. To compute the solution $\boldsymbol{u}$ of the linear system $A\boldsymbol{u} = \boldsymbol{b}$ for a given right hand side $\boldsymbol{b}$, we can proceed as follows:

1. Use the final row to compute $u_{n-1} = b_{n-1}/A_{n-1,n-1}$
2. With the knowledge of $u_{n-1}$, compute $u_{n-2} = \left(b_{n-2} - A_{n-2,n-1}u_{n-1}\right)/A_{n-2,n-2}$
3. With the knowledge of $u_{n-1}$ and $u_{n-2}$, compute $u_{n-3} = \left(b_{n-3} - A_{n-3,n-2}u_{n-2} A_{n-3,n-1}u_{n-1}\right)/A_{n-3,n-3}$
4. Continue in the same way to compute $u_{n-4}, u_{n-5}, \dots, u_0$.

## Tasks
### Pseudocode
Write down the algorithm for solving $A\boldsymbol{u}=\boldsymbol{b}$
### Computational complexity
Show that the solution of an $n\times n$ upper triangular system requires $n^2$ arithmetic operations.
### Tridiagonal matrix
Now let $A$ be a **triangular** matrix with $A_{ij}=0$ for $|i-j|>1$. Here is an example for $n=5$:
$$
\begin{aligned}
A &=  \begin{pmatrix}
1.096 & 0.3391 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} &\textcolor{lightgray}{0} \\
0.0632 & 1.197 & 0.06495 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0}\\
\textcolor{lightgray}{0} & 0.0555 & 1.127 & 0.0008768 & \textcolor{lightgray}{0}\\
\textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 0.2176 & 1.263 & 0.3005\\
\textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 0.3045 & 1.047
\end{pmatrix}
\end{aligned}
$$
* How many operations are required to reduce this to an upper triangular system if the procedure in the lecture is followed, but the zeros are taken into account?
* What is the structure of the resulting upper triangular system?
* How many operations are required to solve this system?
* What is the computational complexity for solving the linear system $A\boldsymbol{u}=\boldsymbol{b}$ is $A$ is a triangular matrix?

## Practicalities
* Write down your solution and save it as a `.pdf` file. You do not have to typeset your solution and can also submit scans of handwritten workings, provided they are legible.
* Make sure you explain your thinking.
* Upload the file to the submission point on moodle


# Exercise 5: PETSc sparse matrices and linear solvers

#### Set: week 9
#### Due: end of week 10

## PETSc matrices
Create two $3\times 3$ sparse PETSc matrices $A$, $B$.

By using suitable functions (see [`petsc4py.Mat` documentation](https://petsc.org/release/petsc4py/reference/petsc4py.PETSc.Mat.html)), compute

* $AB$ 
* $AB^\top$
* $A+B$

##  Gauss Seidel iteration
Instead of $P=D$, we could also use the lower triangular part of $A$ and set $P=D+L$ where
$$
L = \begin{cases}
A_{ij} & \text{if $i<j$} \\
0 & \text{otherwise}
\end{cases}
$$
Convince yourself that for a given vector $\boldsymbol{r}$ the equation $(D+L)z=r$ can be solved row-by-row, i.e. by computing first $\boldsymbol{z}_0 = \boldsymbol{r}_0/A_{00}$, then computing $\boldsymbol{z}_1 = (\boldsymbol{r}_1 - A_{10}\boldsymbol{z}_0)/A_{11}$, $\boldsymbol{z}_2=(\boldsymbol{r}_2 - A_{20}\boldsymbol{z}_0 - A_{21}\boldsymbol{z}_1)/A_{22}$ and so on. The corresponding preconditioner is also known as the successive overrelaxation (SOR) method. It can be chosen by setting `-pc_type sor`. Run the code with this preconditioner - how does the number of iterations change?

## PETSc solver options
For this exercise we consider the $n\times n$ matrix $A$ which is of the following form

$$
A_{ij} = \begin{cases}
2 + h & \text{if $i=j$}\\
-1 & \text{if $j=i\pm 1$ } \\
-1 & \text{if ($i=0$ and $j=n-1$) or ($i=n-1$ and $j=0$)}
\end{cases}
$$

where $h=1/n$. In other words, the entries on the main diagonal are $2+h$ while the entries on the first two sub-diagonals and in the upper right and lower right corner are $-1$.

An example for $n=8$ is shown here:
$$
A = \begin{pmatrix}
  2+h &  -1 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} &  -1\\
 -1 &   2+h &  -1 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0}\\
\textcolor{lightgray}{0} &  -1 &   2+h &  -1 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0}\\
\textcolor{lightgray}{0} & \textcolor{lightgray}{0} &  -1 &   2+h &  -1 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0}\\
\textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} &  -1 &   2+h &  -1 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0}\\
\textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} &  -1 &   2+h &  -1 & \textcolor{lightgray}{0}\\
\textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} &  -1 &   2+h &  -1\\
 -1 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} &  -1 &   2+h
\end{pmatrix}
$$

In this exercise we will use different PETSc solvers and preconditioners to solve the linear system

$$
A \boldsymbol{u} = \boldsymbol{b}
$$

for a given right hand side $\boldsymbol{b}$.

1. Use `PETSc.Mat().createAIJWithArrays()` to create the matrix $A$ in PETSc CSR format
2. Use `PETSc.Vec().createWithArray()` to create a right-hand side vector $\boldsymbol{b}$ which contains random values. You can use the following code to create a numpy array of length $n$ with normally distributed random values:
```Python
rng = np.random.default_rng(seed=1241773)
array = rng.normal(size=n)
```
3. Create a PETSc `KSP` object, which can be configured with PETSc options passed from the command line.
4. Solve the system for different problem sizes $n=32,64,128,256,512$ to a (relative) tolerance of $10^{-9}$ on the preconditioned residual. Investigate the number of iterations and runtime.

Use the following solvers

* Conjugate Gradient (CG)
* Generalised minimal residual (GMRES)
* Richardson iteration

and preconditioners

* Jacobi
* Successive overrelaxation (SOR)
* Algebraic multigrid (`gamg`)
* Incomplete LU factorisation (ILU)

Do all solver/preconditioner combinations work as expected?

Visualise the results.

### Practicalities
* You can extract the number of solver iterations with `ksp.getIterationNumber()`
* If you pass `-ksp_converged_reason `, PETSc will inform you whether the solver has converged.
* You can measure the time spent in the solve by using the `meausure_time` decorator from `fem.utilities`:
```Python
from fem.utilities import measure_time

with measure_time("solve"):
    # call solver here

```
