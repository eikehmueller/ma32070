import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

x = np.array([0.0, 0.2, 0.6, 1.0])
n_points = len(x)
w = np.eye(n_points, n_points)

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
ax = plt.gca()
ax.set_aspect(0.25)
X = np.arange(0, 1, 1e-3)
plt.plot(X, 0 * X, linewidth=1, color="black", linestyle="--")
for j in range(n_points):
    p = sp.interpolate.lagrange(x, w[j, :])
    Y = p(X)
    plt.plot(X, Y, linewidth=2, color=colors[j])
plt.plot(
    x,
    np.zeros(n_points),
    marker="o",
    color="black",
    markerfacecolor="white",
    markersize=8,
    markeredgewidth=2,
    linewidth=0,
)
for j in range(n_points):
    plt.plot(
        [x[j]],
        [1],
        marker="o",
        markerfacecolor="white",
        markersize=6,
        linewidth=0,
        markeredgewidth=2,
        color=colors[j],
    )
plt.savefig("lagrange_polynomials_1d.png", bbox_inches="tight", dpi=300)
