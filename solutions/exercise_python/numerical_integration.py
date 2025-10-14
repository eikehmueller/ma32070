"""Collection of classes for numerical integration"""

import numpy as np
from abc import ABC, abstractmethod

__all__ = ["MidpointRule", "SimpsonsRule"]


class NumericalIntegration(ABC):
    """Numerically integrate a real-valued function in the interval [a,b]

    Numerically approximate the integral

        int_a^b f(x) dx

    by sub-dividing [a,b] into n subintervals
    """

    def __init__(self, interval, n, order=None):
        """Initialise instance

        :arg interval: interval [a,b]
        :arg n: number of subintervals
        :arg order: order of the integrator
        """

        self._a = interval[0]
        self._b = interval[1]
        self._n = n
        self.order = -1 if order is None else order

    def evaluate(self, f):
        """Numerically approximate the integral int_a^b f(x) dx

        For this, loop over the subintervals [x_j,x_{j+1}], approximate
        the integral in each subinterval and sum the result

        :arg f: function to integrate
        """
        h = (self._b - self._a) / self._n
        return float(
            np.sum(
                (
                    self._integrate(f, self._a + j * h, self._a + (j + 1) * h)
                    for j in range(self._n)
                )
            )
        )

    @abstractmethod
    def _integrate(self, f, x_m, x_p):
        """Approximate int_{x_-}^{x_+} f(x)dx

        :arg f: function to integrate
        :arg x_m: lower bound
        :arg x_+: upper bound
        """


class MidpointRule(NumericalIntegration):
    """Numerical integration with the midpoint rule"""

    def __init__(self, interval, n):
        """Initialise instance

        :arg interval: interval [a,b]
        :arg n: number of subintervals
        """
        super().__init__(interval, n, order=2)

    def _integrate(self, f, x_m, x_p):
        """Approximate int_{x_-}^{x_+} f(x)dx by midpoint rule

        (x_+ - x_-) * f((x_+ + x_-)/2)

        :arg f: function to integrate
        :arg x_m: lower bound x_-
        :arg x_p: upper bound x_+
        """
        return (x_p - x_m) * f((x_m + x_p) / 2)


class SimpsonsRule(NumericalIntegration):
    """Numerical integration with Simpson's rule"""

    def __init__(self, interval, n):
        """Initialise instance

        :arg interval: interval [a,b]
        :arg n: number of subintervals
        """
        super().__init__(interval, n, order=4)

    def _integrate(self, f, x_m, x_p):
        """Approximate int_{x_-}^{x_+} f(x)dx by Simpson's rule

        (x_+ - x_-)/6 * ( f(x_-) + 4*f((x_+ + x_-)/2) + f(x_+) )

        :arg f: function to integrate
        :arg x_m: lower bound x_-
        :arg x_p: upper bound x_+
        """
        return (x_p - x_m) / 6 * (f(x_m) + 4 * f((x_m + x_p) / 2) + f(x_p))
