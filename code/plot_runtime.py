import numpy as np
from matplotlib import pyplot as plt


def plot_runtime_error_dense(timing, ndof, error, extension=".pdf"):
    """Plot the runtime and error for dense matrices

    :arg timing: Dictionary with timing data for dense matrices
    :arg ndof: Dictionary with degrees of freedom for each problem size
    :arg error: Dictionary with error values for each problem size
    """
    plt.clf()
    fig, (ax_timing, ax_error) = plt.subplots(1, 2, figsize=(12, 5))
    for label, t in timing.items():
        Y = np.asarray(list(t.values()))
        X = np.asarray(list(ndof.values())[: len(Y)])

        ax_timing.plot(X, Y, linewidth=2, markersize=6, marker="o", label=label)
    ax_timing.plot(
        [1e2, 2e3],
        [0.3e-1, 0.6e0],
        linewidth=2,
        linestyle="--",
        color="black",
        label=r"$\propto n_{\text{dof}}$",
    )
    ax_timing.plot(
        [0.5e3, 4e3],
        [2.0e-4, 8**3 * 2.0e-4],
        linewidth=2,
        linestyle=":",
        color="black",
        label=r"$\propto n_{\text{dof}}^3$",
    )
    ax_timing.set_yscale("log")
    ax_timing.set_xscale("log")
    ax_timing.set_xlabel(r"problem size $n_{\text{dof}}$")
    ax_timing.set_ylabel("time [s]")
    ax_timing.legend(loc="lower right")

    Y = np.asarray(list(error.values()))
    X = np.asarray(list(ndof.values())[: len(Y)])
    ax_error.plot(X, Y, linewidth=2, marker="o", markersize=6)
    ax_error.set_xscale("log")
    ax_error.set_yscale("log")
    ax_error.set_xlabel(r"problem size $n_{\text{dof}}$")
    ax_error.set_ylabel(r"$L_2$ error $\|u^{(h)}-u_{\text{exact}}\|_{L_2(\Omega)}$")
    ax_error.plot(
        [2e2, 5e3],
        [1e-1, 1 / 25 * 1e-1],
        linewidth=2,
        linestyle="--",
        color="black",
        label=r"$\propto n_{\text{dof}}^{-1}\propto h^2$",
    )
    ax_error.legend(loc="upper right")
    plt.savefig("runtime_dense" + extension, bbox_inches="tight")


def plot_niter(niter, ndof, extension=".pdf"):
    """Plot the runtime and error for dense matrices

    :arg niter: Dictionary with number of iterations for each KSP/PC and problem size
    :arg ndof: Dictionary with degrees of freedom for each problem size
    :arg error: Dictionary with error values for each problem size
    """
    plt.clf()
    ndof = list(ndof.values())
    s = r"| $\boldsymbol{n_{\text{dof}}}$ | "
    for label, _ in niter.items():
        ksp, pc = label.split("_")
        s += f"{ksp} + {pc} |"
    print(s)
    print("| " + ((len(niter) + 1) * ":----: |"))
    for j in range(len(ndof)):
        s = "| "
        s += f"{ndof[j]} |"
        for _, _niter in niter.items():
            try:
                _nits = list(_niter.values())[j]
            except:
                _nits = "---"
            s += f" {_nits} |"
        print(s)


def plot_runtime_sparse(timing_sparse, timing_dense, ndof, extension=".pdf"):
    """Plot the runtime and error for sparse matrices

    :arg timing_sparse: Dictionary with timing data for sparse matrices
    :arg timing_dense: Dictionary with timing data for dense matrices
    :arg error: Dictionary with error values for each problem size
    """
    plt.clf()
    ax = plt.gca()
    ax.set_aspect(0.5)
    for solver_label, t_dict in timing_sparse.items():
        ksp, pc = solver_label.split("_")
        label = f"{ksp} + {pc}"
        t_sparse = t_dict["solve()"]
        Y = np.asarray(list(t_sparse.values()))
        X = np.asarray(list(ndof.values())[: len(Y)])
        plt.plot(X, Y, linewidth=2, markersize=6, marker="o", label=label)

    t_dense = timing_dense["solve()"]
    Y = np.asarray(list(t_dense.values()))
    X = np.asarray(list(ndof.values())[: len(Y)])
    plt.plot(X, Y, linewidth=2, markersize=6, marker="o", label="Gaussian elimination")
    plt.plot(
        [0.5e4, 8e4],
        [2.0e-3, 16 * 2.0e-3],
        linewidth=2,
        linestyle="--",
        color="black",
        label=r"$\propto n_{\text{dof}}$",
    )
    plt.plot(
        [2e3, 16e3],
        [1.0e-2, 8**3 * 1.0e-2],
        linewidth=2,
        linestyle=":",
        color="black",
        label=r"$\propto n_{\text{dof}}^3$",
    )
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel(r"problem size $n_{\text{dof}}$")
    ax.set_ylabel("time [s]")
    ax.legend(loc="upper left")

    plt.savefig("runtime_sparse" + extension, bbox_inches="tight")


def plot_time_per_iteration(timing_sparse, niter, ndof, extension=".pdf"):
    """Plot the runtime and error for sparse matrices

    :arg timing_sparse: Dictionary with timing data for sparse matrices
    :arg niter: Dictionary with number of iterations
    :arg error: Dictionary with error values for each problem size
    """
    plt.clf()
    ax = plt.gca()
    ax.set_aspect(0.75)
    for solver_label, t_dict in timing_sparse.items():
        ksp, pc = solver_label.split("_")
        label = f"{ksp} + {pc}"
        t_sparse = t_dict["solve()"]
        Y = np.asarray(list(t_sparse.values())) / np.asarray(
            list(niter[solver_label].values())
        )
        X = np.asarray(list(ndof.values())[: len(Y)])
        plt.plot(X, Y, linewidth=2, markersize=6, marker="o", label=label)
    plt.plot(
        [0.5e4, 8e4],
        [2.0e-4, 16 * 2.0e-4],
        linewidth=2,
        linestyle="--",
        color="black",
        label=r"$\propto n_{\text{dof}}$",
    )
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel(r"problem size $n_{\text{dof}}$")
    ax.set_ylabel("time per iteration [s]")
    ax.legend(loc="upper left")

    plt.savefig("time_per_iteration_sparse" + extension, bbox_inches="tight")


# dense numpy matrices, solve with np.linalg.solve()
timing_dense = {
    "assemble_rhs()": {
        2: 1.29e-03,
        3: 7.99e-03,
        4: 1.57e-02,
        5: 6.43e-02,
        6: 2.57e-01,
        7: 1.13e00,
    },
    "assemble_lhs()": {
        2: 1.40e-03,
        3: 6.76e-03,
        4: 3.01e-02,
        5: 9.70e-02,
        6: 5.00e-01,
        7: 2.19e00,
    },
    "solve()": {
        2: 1.68e-05,
        3: 7.67e-05,
        4: 5.62e-04,
        5: 2.16e-02,
        6: 7.78e-01,
        7: 4.84e01,
    },
}

# PETSc solver -ksp_rtol 1E-9
timing_sparse = {
    "richardson_jacobi": {
        "assemble_rhs()": {
            2: 1.74e-03,
            3: 3.93e-03,
            4: 1.58e-02,
            5: 5.78e-02,
            6: 2.25e-01,
            7: 9.25e-01,
            8: 3.72e00,
            9: 1.47e01,
        },
        "assemble_lhs()": {
            2: 5.91e-03,
            3: 8.75e-03,
            4: 2.13e-02,
            5: 8.26e-02,
            6: 3.27e-01,
            7: 1.31e00,
            8: 5.44e00,
        },
        "solve()": {
            2: 3.66e-03,
            3: 4.93e-03,
            4: 1.81e-02,
            5: 2.21e-01,
            6: 7.70e-01,
            7: 3.34e00,
            8: 2.83e01,
        },
    },
    "cg_jacobi": {
        "assemble_rhs()": {
            2: 2.51e-03,
            3: 4.39e-03,
            4: 1.63e-02,
            5: 5.66e-02,
            6: 2.25e-01,
            7: 9.07e-01,
            8: 3.72e00,
            9: 1.45e01,
        },
        "assemble_lhs()": {
            2: 5.39e-03,
            3: 9.86e-03,
            4: 2.50e-02,
            5: 8.54e-02,
            6: 3.25e-01,
            7: 1.29e00,
            8: 5.36e00,
            9: 2.11e01,
        },
        "solve()": {
            2: 2.53e-03,
            3: 3.27e-03,
            4: 1.79e-03,
            5: 2.61e-03,
            6: 5.88e-03,
            7: 2.18e-02,
            8: 1.36e-01,
            9: 1.02e00,
        },
    },
    "cg_hypre": {
        "assemble_rhs()": {
            2: 1.04e-03,
            3: 3.65e-03,
            4: 1.64e-02,
            5: 5.73e-02,
            6: 2.32e-01,
            7: 9.07e-01,
            8: 3.93e00,
            9: 1.44e01,
        },
        "assemble_lhs()": {
            2: 1.46e-03,
            3: 5.29e-03,
            4: 2.33e-02,
            5: 8.18e-02,
            6: 3.30e-01,
            7: 1.30e00,
            8: 5.36e00,
            9: 2.10e01,
        },
        "solve()": {
            2: 5.13e-03,
            3: 3.17e-04,
            4: 6.00e-04,
            5: 1.57e-03,
            6: 4.67e-03,
            7: 1.92e-02,
            8: 7.08e-02,
            9: 3.15e-01,
        },
    },
}

niter = {
    "richardson_jacobi": {
        2: 3101,
        3: 7959,
        4: 14914,
        5: 32137,
        6: 33888,
        7: 37030,
        8: 129815,
    },
    "cg_jacobi": {
        2: 11,
        3: 23,
        4: 49,
        5: 97,
        6: 187,
        7: 232,
        8: 450,
        9: 876,
    },
    "cg_hypre": {
        2: 3,
        3: 5,
        4: 5,
        5: 6,
        6: 6,
        7: 6,
        8: 6,
        9: 6,
    },
}

error = {
    2: 0.45333811330883717,
    3: 0.1591907906835339,
    4: 0.04709581162800813,
    5: 0.012327845475626148,
    6: 0.0031194424825538424,
    7: 0.0007822955382497838,
}

ndof = {
    2: 25,
    3: 81,
    4: 289,
    5: 1089,
    6: 4225,
    7: 16641,
    8: 66049,
    9: 263169,
}


if __name__ == "__main__":
    extension = ".svg"
    plot_runtime_error_dense(timing_dense, ndof, error, extension=extension)
    plot_niter(niter, ndof, extension=extension)
    plot_runtime_sparse(timing_sparse, timing_dense, ndof, extension=extension)
    plot_time_per_iteration(timing_sparse, niter, ndof, extension=extension)
