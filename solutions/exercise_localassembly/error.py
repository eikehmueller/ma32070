import numpy as np
from fem.quadrature import (
    GaussLegendreQuadratureReferenceTriangle,
)

__all__ = ["error_norm"]


def error_norm(u_numerical, u_exact, element, n_q):
    """Compute L2 norm of the error

    :arg u_numerical: dof-vector of function on reference triangle
    :arg u_exact: exact function which can be evaluated on reference triangle
    :arg element: finite element
    :arg n_q: number of quadrature points
    """
    quad = GaussLegendreQuadratureReferenceTriangle(n_q)
    zeta_q = np.asarray(quad.nodes)
    w_q = quad.weights
    phi = element.tabulate(zeta_q)
    u_q = u_exact(zeta_q.T)
    e_q = u_q - phi @ u_numerical
    nrm2 = np.sum(w_q * e_q**2)
    return np.sqrt(nrm2)
