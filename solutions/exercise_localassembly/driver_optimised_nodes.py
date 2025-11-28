from matplotlib import pyplot as plt
import numpy as np
from firedrake import *


def firedrake_nodes(degree):
    """Return array with Lagrange element nodes used by Firedrake

    :arg degree: polynomial degree
    """

    mesh = UnitTriangleMesh()
    x, y = SpatialCoordinate(mesh)
    V = FunctionSpace(mesh, "CG", degree)
    u_x = Function(V).interpolate(x)
    u_y = Function(V).interpolate(y)

    coordinates = np.asarray([u_x.dat.data, u_y.dat.data]).T
    ndof = V.dof_count
    n_interior = (degree - 2) * (degree - 1) // 2
    indirection_map = (
        [ndof - 3, ndof - 2, ndof - 1]
        + list(range(n_interior + (degree - 1), n_interior + 2 * (degree - 1)))
        + list(range(n_interior + 2 * (degree - 1), n_interior + 3 * (degree - 1)))
        + list(range(n_interior, n_interior + degree - 1))
        + list(range(n_interior))
    )
    return coordinates[indirection_map, :]


def visualise_firedrake_nodes(degree):
    nodes = firedrake_nodes(degree)
    plt.clf()
    for j, X in enumerate(nodes):
        plt.plot([X[0]], [X[1]], linewidth=0, marker="o", markersize=4, color="blue")
        plt.annotate(str(j), X)
    plt.savefig("nodes.pdf", bbox_inches="tight")


"""Main program"""

import functools
import numpy as np


from fem.polynomialelement import PolynomialElement
from fem.utilities import plot_solution

from assembly import (
    assemble_rhs,
    assemble_lhs,
)

from error import error_norm


def u_exact(x, sigma, x0):
    """Analytical solution

    :arg x: point at which the function is evaluated
    :arg sigma: width of peak
    :arg x0: location of peak"""
    return np.exp(-1 / (2 * sigma**2) * ((x[0] - x0[0]) ** 2 + (x[1] - x0[1]) ** 2))


def f(x, kappa, omega, sigma, x0):
    """function to interpolate

    :arg x: point at which the function is evaluated
    :arg kappa: coefficient of diffusion term
    :arg omega: coefficient of zero-order term
    :arg sigma: width of peak
    :arg x0: location of peak
    """
    x_sq = (x[0] - x0[0]) ** 2 + (x[1] - x0[1]) ** 2
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
    if np.all(x[1]) < 1e-12:
        # facet F_1
        n_dot_x = -(x[1] - x0[1])
    elif np.all(x[0]) < 1e-12:
        # facet F_2
        n_dot_x = -(x[0] - x0[0])
    else:
        # facet F_0
        n_dot_x = (x[0] - x0[0] + x[1] - x0[1]) / np.sqrt(2)
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
degree = 10
# Quadrature parameter
n_q = degree + 1

element = PolynomialElement(degree)
# monkey-patch class
element._nodal_points = firedrake_nodes(degree)
vandermonde_matrix = element._vandermonde_matrix(element._nodal_points)
element._coefficients = np.linalg.inv(vandermonde_matrix)

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
nrm = error_norm(
    u_numerical, functools.partial(u_exact, sigma=sigma, x0=x0), element, n_q
)
print(np.linalg.cond(stiffness_matrix))

print(f"{degree}: {nrm:10.4e},")

plot_solution(
    u_numerical,
    functools.partial(u_exact, sigma=sigma, x0=x0),
    element,
    "triangle_solution.pdf",
)
