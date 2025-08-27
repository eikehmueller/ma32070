import numpy as np
from matplotlib import pyplot as plt


def t_matmat(n, p):
    t_flop = 1.77e-11
    t_mem = 2.74e-10
    t_lat = 1.34e-6
    t_word = 2.7e-9
    return (2 * n**3 * t_flop + 3 * n**2 * t_mem) / p + n**2 * t_word + p * t_lat


P = 2 ** np.arange(10)

plt.clf()
fig, axs = plt.subplots(1, 3, figsize=(4, 3))
fig.subplots_adjust(left=0.2, right=3.4)
marker = ["o", "s", "v", "^"]
for j, n in enumerate((256, 1024, 4096, 16384)):
    T = t_matmat(n, P)
    axs[0].plot(
        P,
        T,
        linewidth=2,
        marker=marker[j],
        markersize=6,
        markeredgewidth=2,
        markerfacecolor="white",
    )
    axs[1].plot(
        P,
        t_matmat(n, 1) / T,
        linewidth=2,
        marker=marker[j],
        markersize=6,
        markeredgewidth=2,
        markerfacecolor="white",
        label=f"$n={n}$",
    )
    axs[2].plot(
        P,
        t_matmat(n, 1) / (T * P),
        linewidth=2,
        marker=marker[j],
        markersize=6,
        markeredgewidth=2,
        markerfacecolor="white",
    )
axs[1].plot(P, P, linewidth=2, color="black", linestyle="--")
for ax in axs:
    ax.set_xscale("log")
    ax.set_xlabel("number of processors")
axs[0].set_yscale("log")
axs[0].set_ylabel("runtime $T_p$ [s]")
axs[1].set_ylabel("speedup $T_1/T_p$")
axs[1].set_yscale("log")
axs[1].legend(loc="upper left")
axs[2].set_ylabel(r"parallel efficiency $T_1/(p\cdot T_p)$")
axs[2].set_yticks([0, 0.25, 0.5, 0.75, 1])
axs[2].set_yticklabels([r"0%", r"25%", r"50%", r"75%", r"100%"])
plt.savefig("parallel_scaling.svg", bbox_inches="tight")
