"""Solve finite element model problem and compute global L2 error"""

import numpy as np

from fem.utilitymeshes import rectangle_mesh
from fem.linearelement import LinearElement
from fem.utilities import measure_time
from fem.functionspace import FunctionSpace
from fem.function import Function, CoFunction
from fem.assembly import assemble_rhs, assemble_lhs
from fem.quadrature import GaussLegendreQuadratureReferenceTriangle
from error import error_norm


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
# element = PolynomialElement(3)

# Mesh
mesh = rectangle_mesh(Lx=1, Ly=1, nref=nref)
# Function space
fs = FunctionSpace(mesh, element)
print(f"nref = {nref}")
print(f"grid spacing h = {np.sqrt(2)/2**nref}")
print(f"number of unknowns = {fs.ndof}")
# Quadrature rule
quad = GaussLegendreQuadratureReferenceTriangle(3)

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

error_nrm = error_norm(u_h, u_exact, quad)

print()
print(f"error = {error_nrm}")
