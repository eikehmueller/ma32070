import numpy as np

from cubicelement import CubicElement


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


def test_tabulate_dofs_selected():
    """Check that dof-tabulation of exp(-x)*(2+sin(x)) is correct"""
    element = CubicElement()
    expected = [
        2.00000000,
        0.73575888,
        2.84147098,
        1.19482160,
        1.87614395,
        2.61836980,
        2.32719470,
        1.43306262,
        1.02683424,
        1.66750787,
    ]
    tabulated = element.tabulate_dofs(lambda x: np.exp(-x[0]) * (2 + np.sin(x[1])))
    assert np.allclose(expected, tabulated, rtol=1.0e-8)


def test_tabulate_single_point():
    """Check that tabulation of all dofs at a single point is correct"""
    element = CubicElement()
    zeta = [0.18, 0.43]
    expected = [
        -0.0275145,
        0.060444,
        -0.0442685,
        -0.160218,
        0.101007,
        0.2188485,
        0.1282905,
        0.053703,
        -0.145314,
        0.815022,
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
            -0.0275145,
            0.060444,
            -0.0442685,
            -0.160218,
            0.101007,
            0.2188485,
            0.1282905,
            0.053703,
            -0.145314,
            0.815022,
        ],
        [
            0.0494935,
            0.066816,
            0.0532245,
            0.789264,
            -0.251748,
            -0.0244755,
            -0.0522585,
            -0.179172,
            0.263088,
            0.285768,
        ],
        [
            0.0213005,
            -0.032,
            0.0116095,
            0.1116,
            -0.03906,
            -0.0283185,
            -0.0525915,
            -0.06786,
            0.1044,
            0.97092,
        ],
    ]
    assert np.allclose(tabulated, expected, rtol=1.0e-12)


def test_tabulate_gradient_single_point():
    """Check that tabulation of gradient of all dofs at a single point is correct"""
    element = CubicElement()
    zeta = [0.18, 0.43]
    tabulated = element.tabulate_gradient(zeta)
    expected = [
        [0.45665, 0.45665],
        [-0.1826, 0.0],
        [0.0, -0.37385],
        [0.1548, -0.3726],
        [0.56115, 1.2798],
        [-0.56115, 2.21175],
        [-2.5929, -2.29455],
        [-0.78705, -1.0854],
        [0.513, 0.3726],
        [2.4381, -0.1944],
    ]
    assert np.allclose(tabulated, expected, rtol=1.0e-12)


def test_tabulate_gradient_multiple_points():
    """Check that tabulation of gradient of all dofs at multiple points is correct"""
    element = CubicElement()
    zeta = [[0.18, 0.43], [0.72, 0.21], [0.4, 0.31]]
    tabulated = element.tabulate_gradient(zeta)
    expected = [
        [
            [0.45665, 0.45665],
            [-0.1826, 0.0],
            [0.0, -0.37385],
            [0.1548, -0.3726],
            [0.56115, 1.2798],
            [-0.56115, 2.21175],
            [-2.5929, -2.29455],
            [-0.78705, -1.0854],
            [0.513, 0.3726],
            [2.4381, -0.1944],
        ],
        [
            [-0.43615, -0.43615],
            [1.5184, 0.0],
            [-0.0, -0.29465],
            [3.1374, 3.7584],
            [-0.34965, 0.8424],
            [0.34965, 0.43155],
            [0.5481, 0.29925],
            [1.63035, 1.8792],
            [-2.7126, -3.7584],
            [-3.6855, -2.7216],
        ],
        [
            [0.47465, 0.47465],
            [-0.44, 0.0],
            [0.0, -0.49265],
            [1.953, 0.36],
            [-0.09765, 1.548],
            [0.09765, 1.21995],
            [-1.0323, -1.20195],
            [-1.50165, -1.332],
            [1.467, -0.36],
            [-0.9207, -0.216],
        ],
    ]
    assert np.allclose(tabulated, expected, rtol=1.0e-12)


def test_tabulate_nodal_points():
    """Check that tabulation at nodal points is correct"""
    element = CubicElement()
    xi = [
        [0, 0],
        [1, 0],
        [0, 1],
        [2 / 3, 1 / 3],
        [1 / 3, 2 / 3],
        [0, 2 / 3],
        [0, 1 / 3],
        [1 / 3, 0],
        [2 / 3, 0],
        [1 / 3, 1 / 3],
    ]
    tabulated = element.tabulate(xi)
    assert np.allclose(tabulated, np.eye(10), rtol=1e-12)


def test_tabulate_dofs():
    """Check that dof-tabulation of basis functions is correct"""
    basis_functions = [
        lambda x: 1
        - 11 / 2 * (x[0] + x[1])
        + 9 * (x[0] + x[1]) ** 2
        - 9 / 2 * (x[0] + x[1]) ** 3,
        lambda x: x[0] * (1 - 9 / 2 * x[0] + 9 / 2 * x[0] ** 2),
        lambda x: x[1] * (1 - 9 / 2 * x[1] + 9 / 2 * x[1] ** 2),
        lambda x: -9 / 2 * x[0] * x[1] * (1 - 3 * x[0]),
        lambda x: -9 / 2 * x[0] * x[1] * (1 - 3 * x[1]),
        lambda x: -9 / 2 * x[1] * (1 - x[0] - 4 * x[1] + 3 * x[1] * (x[0] + x[1])),
        lambda x: 9 / 2 * x[1] * (2 - 5 * (x[0] + x[1]) + 3 * (x[0] + x[1]) ** 2),
        lambda x: 9 / 2 * x[0] * (2 - 5 * (x[0] + x[1]) + 3 * (x[0] + x[1]) ** 2),
        lambda x: -9 / 2 * x[0] * (1 - x[1] - 4 * x[0] + 3 * x[0] * (x[0] + x[1])),
        lambda x: 27 * x[0] * x[1] * (1 - x[0] - x[1]),
    ]
    element = CubicElement()
    tabulated = [element.tabulate_dofs(phi) for phi in basis_functions]
    assert np.allclose(tabulated, np.eye(10), rtol=1e-12)


def test_inverse_dofmap():
    """Check that results of inverse_dofmap() are correct"""
    element = CubicElement()
    predicted = {ell: element.inverse_dofmap(ell) for ell in range(10)}
    expected = {
        0: ("vertex", 0, 0),
        1: ("vertex", 1, 0),
        2: ("vertex", 2, 0),
        3: ("facet", 0, 0),
        4: ("facet", 0, 1),
        5: ("facet", 1, 0),
        6: ("facet", 1, 1),
        7: ("facet", 2, 0),
        8: ("facet", 2, 1),
        9: ("interior", 0, 0),
    }
    assert predicted == expected
