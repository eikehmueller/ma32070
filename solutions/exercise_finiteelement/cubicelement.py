"""Linear finite element"""

import numpy as np
from fem.finiteelement import FiniteElement

__all__ = ["CubicElement"]


class CubicElement(FiniteElement):
    """Cubic finite element basis functions in 2d

    There are 10 basis functions.

    Arrangement of the 10 unknowns on reference triangle:

    V 2
        2
        ! .
        !  .
        5   4
        !    .
        !     . F 0
    F 1 !      .
        !   9   .
        6        3
        !         .
        !          .
        0---7---8---1
    V 0       F 2     V 1


    """

    def __init__(self):
        """Initialise new instance"""
        super().__init__()
        self._nodal_points = np.asarray(
            [
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
        )
        vandermonde_matrix = self._vandermonde_matrix(self._nodal_points)
        # Solve A.C = Id for the coefficient matrix C. The k-th column of C contains
        # the polynomial coefficients for the k-th basis function
        self._coefficients = np.linalg.inv(vandermonde_matrix)

    def _vandermonde_matrix(self, zeta, grad=False):
        """Construct the Vandermonde matrix or its gradient

        If grad=False, compute the Vandermonde matrix V(zeta)

            V_{i,j}(zeta) = x_i^a(j)*y_i^b(j)

        where the row i corresponds to the index of the point zeta_i = (x_i,y_i) and
        the column j to the power (a(j),b(j)) that the point zeta_i is raised to.
        The resulting matrix has the shape (npoints, ndof) where npoints is the number of points
        in zeta.

        If grad=True, compute the gradient grad V(zeta) of the Vandermonde matrix with

            grad V_{i,j,k} = d V_{i,j}(zeta) / dx_k

        The resulting tensor has the shape (npoints, ndof,2).

        :arg zeta: array of shape (npoints, 2) containing the points for which the Vandermonde matrix
                 is to be calculated
        :arg grad: compute gradient?
        """

        npoints = zeta.shape[0]
        if grad:
            mat = np.zeros([npoints, 10, 2])
            mat[:, 1, 0] = 1
            mat[:, 3, 0] = 2 * zeta[..., 0]
            mat[:, 4, 0] = zeta[..., 1]
            mat[:, 6, 0] = 3 * zeta[..., 0] ** 2
            mat[:, 7, 0] = 2 * zeta[..., 0] * zeta[..., 1]
            mat[:, 8, 0] = zeta[..., 1] ** 2

            mat[:, 2, 1] = 1
            mat[:, 4, 1] = zeta[..., 0]
            mat[:, 5, 1] = 2 * zeta[..., 1]
            mat[:, 7, 1] = zeta[..., 0] ** 2
            mat[:, 8, 1] = 2 * zeta[..., 0] * zeta[..., 1]
            mat[:, 9, 1] = 3 * zeta[..., 1] ** 2
        else:
            mat = np.empty([npoints, 10])
            mat[:, 0] = 1
            mat[:, 1] = zeta[..., 0]
            mat[:, 2] = zeta[..., 1]
            mat[:, 3] = zeta[..., 0] ** 2
            mat[:, 4] = zeta[..., 0] * zeta[..., 1]
            mat[:, 5] = zeta[..., 1] ** 2
            mat[:, 6] = zeta[..., 0] ** 3
            mat[:, 7] = zeta[..., 0] ** 2 * zeta[..., 1]
            mat[:, 8] = zeta[..., 0] * zeta[..., 1] ** 2
            mat[:, 9] = zeta[..., 1] ** 3
        return mat

    @property
    def ndof_per_interior(self):
        """Return number of unknowns associated with the interior of the cell"""
        return 1

    @property
    def ndof_per_facet(self):
        """Return number of unknowns associated with each facet"""
        return 2

    @property
    def ndof_per_vertex(self):
        """Return number of unknowns associated with each vertex"""
        return 1

    def tabulate_dofs(self, fhat):
        """Evaluate the dofs on a given function on the reference element

        :arg fhat: function fhat defined for 2d vectors
        """
        return fhat(self._nodal_points.T)

    def tabulate(self, zeta):
        """Evaluate all basis functions at a point inside the reference cell

        Returns a vector of length ndof with the evaluation of all basis functions or a matrix
        of shape (npoints,ndof) if xi contains several points.

        :arg zeta: point zeta=(x,y) at which the basis functions are to be evaluated; can also be a
                 matrix of shape (npoints,2).
        """
        _zeta = np.asarray(zeta)
        mat = np.squeeze(
            self._vandermonde_matrix(
                np.expand_dims(_zeta, axis=list(range(2 - _zeta.ndim))),
                grad=False,
            )
            @ self._coefficients
        )
        return mat

    def tabulate_gradient(self, zeta):
        """Evaluate the gradients of all basis functions at a point inside the reference cell

        Returns an matrix of shape (ndof,2) with the evaluation of the gradients of all
        basis functions. If zeta is a matrix containing several points then the matrix that is
        returned is of shape (npoints,ndof,2)

        :arg zeta: point zeta=(x,y) at which the gradients of the  basis functions are to be evaluated;
                 can also be a matrix of shape (npoints,2).
        """
        _zeta = np.asarray(zeta)
        mat = np.squeeze(
            np.einsum(
                "imk,mj->ijk",
                self._vandermonde_matrix(
                    np.expand_dims(_zeta, axis=list(range(2 - _zeta.ndim))),
                    grad=True,
                ),
                self._coefficients,
            )
        )
        return mat
