import functools
import numpy as np
from numerical_integration import MidpointRule, SimpsonsRule


def f(x, alpha):
    """Function to integrate f(x,alpha) = exp(-alpha*x)"""
    return np.exp(-alpha * x)


alpha = 0.4
exact_result = 1 / alpha * (1 - np.exp(-alpha))

# pedestrian implementation
results = {"MidpointRule": [], "SimpsonsRule": []}
for n in [4, 8, 16, 32]:
    integrator_midpoint = MidpointRule([0, 1], n)
    results["MidpointRule"].append(
        integrator_midpoint.evaluate(functools.partial(f, alpha=alpha))
    )
for n in [4, 8, 16, 32]:
    integrator_simpson = SimpsonsRule([0, 1], n)
    results["SimpsonsRule"].append(
        integrator_simpson.evaluate(functools.partial(f, alpha=alpha))
    )

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
    print(
        f"{integrator}: "
        + ", ".join([f"{abs(x-exact_result):8.4e}" for x in integrals])
    )
