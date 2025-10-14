import functools
import numpy as np
from numerical_integration import MidpointRule, SimpsonsRule


def f(x, alpha):
    """Function to integrate f(x,alpha) = exp(-alpha*x)"""
    return np.exp(-alpha * x)


alpha = 0.4
I_exact = 1 / alpha * (1 - np.exp(-alpha))

# using list- and dictionary comprehensions
results = {
    Integrator.__name__: [
        Integrator([0, 1], n).evaluate(functools.partial(f, alpha=alpha))
        for n in [4, 8, 16, 32]
    ]
    for Integrator in [MidpointRule, SimpsonsRule]
}

# Print out results
for integrator, integrals in results.items():
    print(f"{integrator}: " + ", ".join([f"{abs(I-I_exact):8.4e}" for I in integrals]))

# Print out empirical comvergence rate
for integrator, integrals in results.items():
    print(f"{integrator}: ")
    mu = [
        np.log2(
            (integrals[j] - integrals[j + 1]) / (integrals[j + 1] - integrals[j + 2])
        )
        for j in range(len(integrals) - 2)
    ]
    print(
        "  mu = ",
        ", ".join([f"{z:8.4f}" for z in mu]),
    )
