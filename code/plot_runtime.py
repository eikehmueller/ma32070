import numpy as np
from matplotlib import pyplot as plt

timing = {
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
}

plt.clf()
fig, (ax_timing, ax_error) = plt.subplots(1, 2, figsize=(12, 5))
for label, t in timing.items():
    X = np.asarray(list(ndof.values()))
    Y = np.asarray(list(t.values()))
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

X = np.asarray(list(ndof.values()))
Y = np.asarray(list(error.values()))
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
plt.savefig("runtime.svg", bbox_inches="tight")
