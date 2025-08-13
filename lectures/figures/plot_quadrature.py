import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


npoints = 2
zeta_x, weights_x = np.polynomial.legendre.leggauss(npoints + 1)
zeta_y, weights_y = np.polynomial.legendre.leggauss(npoints)
zeta = []
zeta_Q = []
for x, w_x in zip(zeta_x, weights_x):
    for y, w_y in zip(zeta_y, weights_y):
        zeta_Q.append([x, y])
        x_sq = 1 / 2 * (1 + x)
        y_sq = 1 / 2 * (1 + y)
        zeta.append([x_sq, (1 - x_sq) * y_sq])
zeta = np.asarray(zeta)
zeta_Q = np.asarray(zeta_Q)


plt.clf()
fig, axs = plt.subplots(1, 2)
mask_Q = Polygon(
    [[-1, -1], [+1, -1], [+1, +1], [-1, +1]],
    color="lightgray",
    linewidth=2,
)

ax = axs[0]
ax.set_aspect("equal")
ax.set_xlim(-1.1, 1.1)
ax.set_ylim(-1.1, 1.1)
for x in zeta_x:
    ax.plot([x, x], [-1, 1], linewidth=2, color="gray", linestyle="--")
for y in zeta_y:
    ax.plot([-1, 1], [y, y], linewidth=2, color="gray", linestyle="--")

p = PatchCollection([mask_Q], color="lightgray")
ax.add_collection(p)
ax.plot(
    zeta_Q[:, 0],
    zeta_Q[:, 1],
    linewidth=0,
    marker="o",
    markersize=4,
    color="red",
)
mask = Polygon(
    [[0, 0], [1, 0], [0, 1]],
    color="lightgray",
    linewidth=2,
)
ax.annotate(r"$S$", (0.3, 0), fontsize=24)

ax = axs[1]
ax.set_aspect("equal")
ax.set_xlim(-0.1, 1.1)
ax.set_ylim(-0.1, 1.1)
for x in zeta_x:
    ax.plot(
        [(1 + x) / 2, (1 + x) / 2],
        [0, (1 - x) / 2],
        linewidth=2,
        color="gray",
        linestyle="--",
    )
for y in zeta_y:
    ax.plot([0, 1], [(1 + y) / 2, 0], linewidth=2, color="gray", linestyle="--")

p = PatchCollection([mask], color="lightgray")
ax.add_collection(p)
ax.plot(
    zeta[:, 0],
    zeta[:, 1],
    linewidth=0,
    marker="o",
    markersize=4,
    color="blue",
)
ax.annotate(r"$\widehat{K}$", (0.6, 0.6), fontsize=24)
plt.savefig("quadrature.svg", bbox_inches="tight")
