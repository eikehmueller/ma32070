import numpy as np
import pytest

from threepointquadrature import ThreePointQuadratureReferenceTriangle


@pytest.mark.parametrize(
    "s, expected",
    [
        [[0, 0], 1 / 2],
        [[1, 0], 1 / 6],
        [[0, 1], 1 / 6],
        [[2, 0], 1 / 12],
        [[1, 1], 1 / 24],
        [[0, 2], 1 / 12],
    ],
)
def test_threepoint_quadrature_monomial(s, expected):
    """Check that three point quadrature is exact for monomials x_0^{s_0} x_1^{s_1}
    :arg s: (s_0,s_1) = powers of x_0,x_1
    :arg expected: exact result s_0! s_1! / (s_0+s_1+2)!
    """
    quadrature = ThreePointQuadratureReferenceTriangle()
    s_numerical = 0
    for w, xi in zip(quadrature.weights, quadrature.nodes):
        s_numerical += w * xi[0] ** s[0] * xi[1] ** s[1]

    tolerance = 1.0e-12
    assert np.allclose(s_numerical, expected, tolerance)
