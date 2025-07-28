"""Three point quadrature on reference triangle"""

import numpy as np
from fem.quadrature import Quadrature

__all__ = [
    "ThreePointQuadratureReferenceTriangle",
]


class ThreePointQuadratureReferenceTriangle(Quadrature):
    """Three point quadrature on reference triangle"""

    def __init__(self):
        """Initialise a new instance"""
        super().__init__()
        self._nodes = np.asarray([[1 / 6, 1 / 6], [2 / 3, 1 / 6], [1 / 6, 2 / 3]])
        self._weights = np.asarray([1 / 6, 1 / 6, 1 / 6])

    @property
    def nodes(self):
        """Return quadrature nodes"""
        return self._nodes

    @property
    def weights(self):
        """Return quadrature weights"""
        return self._weights

    @property
    def degree_of_precision(self):
        """Degree of precision, i.e. highest polynomial degree that
        can be integrated exactly"""
        return 2
