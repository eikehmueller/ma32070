import numpy as np
from fem.quadrature import (
    GaussLegendreQuadratureLineSegment,
    GaussLegendreQuadratureReferenceTriangle,
)

__all__ = [
    "assemble_lhs",
    "assemble_rhs",
]


def assemble_lhs(element, n_q, kappa, omega):
    """Assemble LHS bilinear form into matrix

    :arg element: finite element
    :arg n_q: number of quadrature points
    :arg kappa: coefficient kappa of diffusion term
    :arg omega: coefficient kappa of zero-order term
    """
    quad = GaussLegendreQuadratureReferenceTriangle(n_q)
    # extract quadrature points and weights
    zeta_q = np.asarray(quad.nodes)
    w_q = quad.weights
    # tabulation of basis functions and their gradients at quadrature points
    T = element.tabulate(zeta_q)
    T_grad = element.tabulate_gradient(zeta_q)
    stiffness_matrix = kappa * np.einsum(
        "q,qka,qla->lk", w_q, T_grad, T_grad
    ) + omega * np.einsum(
        "q,qk,ql->lk",
        w_q,
        T,
        T,
    )
    return stiffness_matrix


def assemble_rhs(f, g, element, n_q):
    """Assemble functions into RHS on unit triangle

    :arg f: RHS function, needs to be callable with f(x) where x is a 2-vector
    :arg g: Neumann boundary function, needs to be callable with f(x) where x is a 2-vector
    :arg element: finite element
    :arg n_q: number of quadrature points
    """
    quad = GaussLegendreQuadratureReferenceTriangle(n_q)
    # extract quadrature points and weights
    zeta_q = np.asarray(quad.nodes)
    w_q = quad.weights
    f_q = f(zeta_q.T)
    # tabulation of basis functions at quadrature points
    T = element.tabulate(zeta_q)
    r = np.einsum("q,q,ql->l", w_q, f_q, T)
    for v_a, v_b in [[(1, 0), (0, 1)], [(0, 1), (0, 0)], [(0, 0), (1, 0)]]:
        # quadrature rule on facet
        quad_facet = GaussLegendreQuadratureLineSegment(v_a, v_b, n_q)
        # extract quadrature points and weights
        w_facet_q = quad_facet.weights
        zeta_facet_q = np.asarray(quad_facet.nodes)
        # evaluation of function g at quadrature points
        g_q = g(zeta_facet_q.T)
        # tabulation of basis functions at (facet) quadrature points
        T_facet = element.tabulate(zeta_facet_q)
        r += np.einsum("q,q,ql->l", w_facet_q, g_q, T_facet)
    return r
