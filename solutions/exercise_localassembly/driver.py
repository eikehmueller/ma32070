"""Main program"""

import functools
import numpy as np


from fem.polynomialelement import PolynomialElement
from fem.utilities import plot_solution

from algorithms import (
    assemble_rhs,
    assemble_lhs,
    error_nrm,
)


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


# width of peak
sigma = 0.5
# location of peak
x0 = np.asarray([0.6, 0.25])
# Coefficient of diffusion term
kappa = 0.9
# Coefficient of zero order term
omega = 0.4
# Polynomial degree
degree = 1
# Quadrature parameter
n_q = degree + 1

element = PolynomialElement(degree)

# Assemble stiffness matrix
stiffness_matrix = assemble_lhs(element, n_q, kappa, omega)

# Assemble right hand side
rhs = assemble_rhs(
    functools.partial(f, kappa=kappa, omega=omega, sigma=sigma, x0=x0),
    functools.partial(g, kappa=kappa, sigma=sigma, x0=x0),
    element,
    n_q,
)

# Compute solution vector by solving the linear system
u_numerical = np.linalg.solve(stiffness_matrix, rhs)

# Compute the L2 norm of the error ||u - u_exact||_L2
nrm = error_nrm(
    u_numerical, functools.partial(u_exact, sigma=sigma, x0=x0), element, n_q
)

print(f"{degree}: {nrm:10.4e},")

plot_solution(
    u_numerical,
    functools.partial(u_exact, sigma=sigma, x0=x0),
    element,
    "triangle_solution.pdf",
)
