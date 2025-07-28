import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


def plot_solution(u_numerical, u_exact, element, filename):
    """Plot numerical solution, exact solution and error

    :arg u_numerical: numerical solution vector
    :arg u_exact: exact solution function
    :arg element: finite element
    :arg filename: name of file to save plot to
    """

    h = 0.01
    X = np.arange(0, 1 + h / 2, h)
    Y = np.arange(0, 1 + h / 2, h)
    X, Y = np.meshgrid(X, Y)

    XY = np.asarray([X, Y]).T.reshape([X.shape[0] * X.shape[1], 2])

    T = element.tabulate(XY)
    Z_exact = u_exact(XY).reshape([X.shape[0], X.shape[1]]).T
    Z_numerical = np.dot(T, u_numerical).reshape([X.shape[0], X.shape[1]]).T

    fig, axs = plt.subplots(1, 3)
    cs = axs[0].contourf(X, Y, Z_exact, levels=100, vmin=0, vmax=1)
    cbar = fig.colorbar(cs, shrink=1, location="bottom")
    cbar.ax.tick_params(labelsize=4)
    axs[0].set_title("exact")

    cs = axs[2].contourf(X, Y, Z_numerical, levels=100, vmin=0, vmax=1)
    cbar = fig.colorbar(cs, shrink=1, location="bottom")
    cbar.ax.tick_params(labelsize=4)
    axs[2].set_title("numerical")

    cs = axs[1].contourf(
        X, Y, np.log(np.abs(Z_numerical - Z_exact)), levels=100, norm="linear"
    )
    cbar = fig.colorbar(cs, shrink=1, location="bottom")
    cbar.ax.tick_params(labelsize=4)
    axs[1].set_title("log(error)")

    for ax in axs:
        ax.set_aspect("equal")
        border = 2
        mask = Polygon(
            [[1 + border, -border], [-border, 1 + border], [1 + border, 1 + border]],
            color="white",
            linewidth=4,
        )
        p = PatchCollection([mask], color="white", zorder=2)
        ax.add_collection(p)
    plt.savefig(filename, bbox_inches="tight")
