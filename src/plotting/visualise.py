"""Visualise mesh, finite element and quadrature"""

import numpy as np

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

from fem.polynomialelement import PolynomialElement
from fem.utilitymeshes import RectangleMesh
from fem.quadrature import GaussLegendreQuadratureReferenceTriangle


def visualise_mesh(mesh, filename):
    """Plot connectivity information for a mesh

    :arg mesh: mesh object
    :arg filename: name of file to save output to
    """
    plt.clf()
    _, axs = plt.subplots(2, 2)
    plt.subplots_adjust(top=0.99, bottom=0.01, hspace=0.25, wspace=0.4)
    for j in range(2):
        for k in range(2):
            axs[j, k].set_aspect("equal")
    # vertices
    for cell in range(mesh.ncells):
        p = np.zeros((4, 3))
        for j in range(3):
            p[j, :2] = np.asarray(mesh.vertices[mesh.cell2vertex[cell][j]])
        p[-1, :] = p[0, :]
        axs[0, 0].plot(
            p[:, 0],
            p[:, 1],
            linewidth=2,
            color="lightgray",
        )
    axs[0, 0].plot(
        mesh.vertices[:, 0],
        mesh.vertices[:, 1],
        linewidth=0,
        marker="o",
        markersize=4,
        color="blue",
    )
    for vertex in range(mesh.nvertices):
        axs[0, 0].annotate(
            f"{vertex:3d}",
            (mesh.vertices[vertex, :]),
            fontsize=6,
        )
        axs[0, 0].set_title("global vertex index")
    # facets
    for facet in range(mesh.nfacets):
        p = np.asarray([mesh.vertices[vertex] for vertex in mesh.facet2vertex[facet]])
        m = np.mean(p, axis=0)
        rho = 0.8
        p = rho * p + (1 - rho) * np.expand_dims(m, axis=0)
        axs[0, 1].plot(
            p[:, 0], p[:, 1], linewidth=2, marker="o", markersize=4, color="blue"
        )
        axs[0, 1].annotate(
            f"{facet:3d}",
            (m[0], m[1]),
            fontsize=6,
        )
        axs[0, 1].arrow(
            p[0, 0],
            p[0, 1],
            m[0] - p[0, 0],
            m[1] - p[0, 1],
            width=0,
            linewidth=0,
            head_width=0.05,
            color="blue",
        )
    axs[0, 1].set_title("global facet index")
    # cells
    for cell in range(mesh.ncells):
        p = np.zeros((4, 3))
        for j in range(3):
            p[j, :2] = np.asarray(mesh.vertices[mesh.cell2vertex[cell][j]])
        p[-1, :] = p[0, :]
        m = np.mean(p[:-1, :], axis=0)
        rho = 0.8
        p = rho * p + (1 - rho) * np.expand_dims(m, axis=0)
        axs[1, 0].plot(p[:, 0], p[:, 1], linewidth=2, color="blue")
        axs[1, 0].annotate(
            f"{cell:3d}",
            (m[0], m[1]),
            verticalalignment="center",
            horizontalalignment="center",
            fontsize=6,
        )
    axs[1, 0].set_title("global cell index")
    # local indices
    for cell in range(mesh.ncells):
        p = np.zeros((4, 3))
        p_facet = np.zeros((3, 2))
        for j, facet in enumerate(mesh.cell2facet[cell]):
            p_facet[j, :] = 0.5 * (
                mesh.vertices[mesh.facet2vertex[facet][0]]
                + mesh.vertices[mesh.facet2vertex[facet][1]]
            )
            p[j, :2] = np.asarray(mesh.vertices[mesh.cell2vertex[cell][j]])
        p[-1, :] = p[0, :]
        m = np.mean(p[:-1, :], axis=0)
        rho = 0.85
        p_cell = rho * p + (1 - rho) * np.expand_dims(m, axis=0)
        axs[1, 1].plot(
            p_cell[:, 0],
            p_cell[:, 1],
            linewidth=2,
            color="black",
            markersize=4,
            marker="o",
            markerfacecolor="red",
        )
        axs[1, 1].annotate(
            f"{cell:3d}",
            (m[0], m[1]),
            verticalalignment="center",
            horizontalalignment="center",
            fontsize=6,
        )
        omega = 0.6
        p_vertex = omega * p + (1 - omega) * np.expand_dims(m, axis=0)
        p_facet = omega * p_facet + (1 - omega) * np.expand_dims(m[:2], axis=0)
        for j in range(3):
            axs[1, 1].annotate(
                f"{j:d}",
                (p_vertex[j, 0], p_vertex[j, 1]),
                verticalalignment="center",
                horizontalalignment="center",
                color="red",
                fontsize=6,
            )
            axs[1, 1].annotate(
                f"{j:d}",
                p_facet[j, :],
                verticalalignment="center",
                horizontalalignment="center",
                color="blue",
                fontsize=6,
            )
            axs[1, 1].arrow(
                p_cell[j, 0],
                p_cell[j, 1],
                0.5 * (p_cell[(j + 1) % 3, 0] - p_cell[j, 0]),
                0.5 * (p_cell[(j + 1) % 3, 1] - p_cell[j, 1]),
                width=0,
                linewidth=1,
                head_width=0.025,
                color="black",
            )
    axs[1, 1].set_title("local vertex- and facet-index")
    plt.savefig(filename, bbox_inches="tight")


def balanced_factorisation(n):
    """factorise n = a*b such that a+b is minimal and a<b

    :arg n: number to factorise
    """
    s = dict()
    for a in range(1, n + 1):
        if (n // a) * a == n:
            s[(a, n // a)] = a + n // a
    return min(s, key=s.get)


def visualise_element(element, filename):
    """Visualise the basis functions of a finite element

    :arg element: finite element
    :arg filename: name of file that the visualisation is saved to
    """

    plt.clf()
    ndof = element.ndof
    nrows, ncols = balanced_factorisation(ndof)
    fig, axs = plt.subplots(nrows, ncols, figsize=(10, 5))
    h = 0.002
    X = np.arange(0, 1 + h / 2, h)
    Y = np.arange(0, 1 + h / 2, h)
    X, Y = np.meshgrid(X, Y)

    xi = np.asarray((X.flatten(), Y.flatten())).T
    Z = element.tabulate(xi).reshape((*X.shape, ndof))

    h_c = 0.1
    X_c = np.arange(0, 1 + h_c / 2, h_c)
    Y_c = np.arange(0, 1 + h_c / 2, h_c)
    X_c, Y_c = np.meshgrid(X_c, Y_c)
    xi_c = np.asarray((X_c.flatten(), Y_c.flatten())).T
    gradZ = element.tabulate_gradient(xi_c).reshape((*X_c.shape, ndof, 2))
    idx = np.where(
        np.logical_or(
            np.logical_or(X.flatten()[...] > 1, Y.flatten()[...] > 1),
            X.flatten()[...] + Y.flatten()[...] > 1,
        )
    )

    for j in range(ndof):
        row = j // ncols
        col = j % ncols
        if nrows == 1:
            ax = axs[col]
        else:
            ax = axs[row, col]
        if col > 0:
            ax.set_yticks([])
        if row < nrows - 1:
            ax.set_xticks([])
        ax.set_aspect(1)
        ax.set_xlim(-0.1, 1.1)
        ax.set_ylim(-0.1, 1.1)
        ax.set_title(j)
        Z[..., j].reshape([-1])[idx] = 0
        Z[-1, -1] = 1
        Z[-1, -2] = -0.125
        cs = ax.contourf(
            X, Y, Z[..., j], levels=100, vmin=-0.125, vmax=1.0, cmap="terrain"
        )
        ax.quiver(X_c, Y_c, gradZ[..., j, 0], gradZ[..., j, 1], color="red")
        sigma = 2
        mask = Polygon(
            [[1 + sigma, -sigma], [-sigma, 1 + sigma], [1 + sigma, 1 + sigma]],
            color="white",
            linewidth=4,
        )
        p = PatchCollection([mask], color="white", zorder=2)
        ax.add_collection(p)

        ax.plot(
            element._nodal_points[:, 0],
            element._nodal_points[:, 1],
            linewidth=0,
            markersize=4,
            markerfacecolor="white",
            marker="o",
            color="red",
        )
        ax.plot(
            element._nodal_points[j, 0],
            element._nodal_points[j, 1],
            linewidth=0,
            markersize=4,
            marker="o",
            color="red",
        )
    fig.colorbar(cs, ax=axs[...], shrink=0.6, location="right", extend="both")

    plt.savefig(filename, bbox_inches="tight")


def visualise_quadrature(quad, filename):
    """Visualise quadrature rule

    :arg quad: quadrature rule to visualise
    :arg filename: name of file to write to
    """
    plt.clf()
    mask = Polygon(
        [[0, 0], [1, 0], [0, 1]],
        color="lightgray",
        linewidth=2,
    )
    ax = plt.gca()
    ax.set_aspect("equal")
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.1)
    p = PatchCollection([mask], color="lightgray")
    ax.add_collection(p)
    plt.plot(
        quad.nodes[:, 0],
        quad.nodes[:, 1],
        linewidth=0,
        marker="o",
        markersize=1,
        color="blue",
    )
    plt.savefig(filename, bbox_inches="tight")


def to_latex(matrix):
    """Convert a numpy matrix to a latex string

    :arg matrix: Matrix to convert
    """
    _matrix = np.asarray(matrix)
    nrow, ncol = _matrix.shape
    s = ""
    s += r"\begin{pmatrix}" + "\n"
    for j in range(nrow):
        s += " & ".join([str(x) for x in _matrix[j, :]])
        if j < nrow - 1:
            s += r"\\"
        s += "\n"
    s += r"\end{pmatrix}"
    return s


mesh = RectangleMesh(Lx=1.0, Ly=1.0, nref=1)

print("cell2facet:")
print(to_latex(np.asarray(mesh.cell2facet).T))
print()
print("facet2vertex:")
print(to_latex(np.asarray(mesh.facet2vertex).T))
print()
print("cell2vertex:")
print(to_latex(np.asarray(mesh.cell2vertex).T))
print()

visualise_mesh(mesh, "mesh.svg")

# try:
from fem.polynomialelement import PolynomialElement

element = PolynomialElement(2)
visualise_element(element, "element.png")
# except:
#    print("FAIL")

quad = GaussLegendreQuadratureReferenceTriangle(3)
visualise_quadrature(quad, "quadrature.pdf")
