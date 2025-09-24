import numpy as np
from matplotlib import pyplot as plt

grid_spacing = [
    0.1767766952966369,
    0.08838834764831845,
    0.04419417382415922,
    0.02209708691207961,
    0.011048543456039806,
]

error = {
    "LinearElement": [
        0.15919079068353392,
        0.04709581162800816,
        0.012327845475626155,
        0.0031194424825538407,
        0.0007822955382497918,
    ],
    "CubicElement": [
        0.0019296912427000641,
        0.00011906213963505698,
        7.3308506308042225e-06,
        4.556608092935774e-07,
    ],
}

plt.clf()
ax = plt.gca()
for element, v in error.items():
    X = np.asarray(grid_spacing[: len(v)])
    Y = np.asarray(v)
    plt.plot(
        X,
        Y,
        markersize=6,
        marker="o",
        markerfacecolor="white",
        markeredgewidth=2,
        linewidth=2,
        label=element,
    )

rho = 4
x0 = 2e-2
y0 = 6e-3

plt.plot(
    [x0, rho * x0],
    [y0, rho**2 * y0],
    linewidth=2,
    color="black",
    linestyle="--",
    label=r"$\propto h^{2}$",
)

rho = 4
x0 = 3e-2
y0 = 4e-6

plt.plot(
    [x0, rho * x0],
    [y0, rho**4 * y0],
    linewidth=2,
    color="black",
    linestyle=":",
    label=r"$\propto h^{4}$",
)

ax.set_xlabel("grid spacing $h$")
ax.set_ylabel(r"error norm $\|e_h\|_{L_2(\Omega)}$")
ax.set_xscale("log")
ax.set_yscale("log")
plt.legend(loc="lower right")
plt.savefig("convergence.png", bbox_inches="tight", dpi=300)
