import numpy as np
from matplotlib import pyplot as plt

t_flop = 1.6e-11
t_mem = 3.2e-10


def performance(q):
    return 1 / t_flop * 1 / (1 + t_mem / (q * t_flop))


logQ = np.arange(-4, 12, 0.01)

Q = np.exp(logQ)

plt.clf()
plt.plot(
    Q,
    performance(Q),
    linestyle="--",
    color="black",
    label=r"$R_{\text{flop}} \frac{1}{1 + \frac{1}{q} \frac{t_{\text{mem}}}{t_{\text{flop}}}}$",
)
plt.plot(Q, Q / t_mem, color="black")
plt.plot(Q, Q / Q / t_flop, color="black")
plt.fill_between(Q, 1e14, performance(Q), color="lightgray", alpha=0.5)
ax = plt.gca()
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim(1e8, 5e12)
ax.set_xlim(np.exp(np.min(logQ)), np.exp(np.max(logQ)))
plt.text(
    0.04,
    1.5 / t_flop,
    r"$R_{\text{flop}} = t_{\text{flop}}^{-1}$",
    fontsize=12,
    color="black",
)

plt.text(
    21,
    42 / t_mem,
    r"slope = $BW = t_{\text{mem}}^{-1}$",
    fontsize=12,
    rotation=50,
    color="black",
)
plt.text(
    1.5,
    0.125e9,
    r"bandwidth bound" + "\n" + r"$q\ll t_{\text{mem}}/t_{\text{flop}}$",
    ha="center",
    fontsize=12,
)
plt.text(
    1.5e4,
    0.125e9,
    r"compute bound" + "\n" + r"$q\gg t_{\text{mem}}/t_{\text{flop}}$",
    ha="center",
    fontsize=12,
)

plt.text(
    1.5 * t_mem / t_flop,
    0.2e9,
    r"$\frac{t_{\text{mem}}}{t_{\text{flop}}}$",
    ha="left",
    color="blue",
    fontsize=18,
)

plt.text(0.8e2, 2e10, "good code", color="green", fontsize=12)
plt.text(0.5e2, 3e9, "bad code", color="red", fontsize=12)
plt.plot(
    [1.7e3, 1e4],
    [2.5e10, 0.75 * performance(1e4)],
    color="green",
)
plt.plot(
    [0.6e2, 2],
    [2e10, 0.75 * performance(2)],
    color="green",
)
plt.plot(
    [0.8e3, 1.5e4],
    [3.5e9, 0.1 * performance(1.5e4)],
    color="red",
)

plt.plot(
    [0.4e2, 3],
    [3.5e9, 0.15 * performance(3)],
    color="red",
)


plt.plot(
    [2, 1e4],
    [0.75 * performance(2), 0.75 * performance(1e4)],
    linewidth=0,
    marker="o",
    markeredgecolor="green",
    markeredgewidth=2,
    markerfacecolor="white",
    markersize=10,
)
plt.plot(
    [3, 1.5e4],
    [0.15 * performance(3), 0.1 * performance(1.5e4)],
    color="red",
    linewidth=0,
    marker="o",
    markeredgecolor="red",
    markeredgewidth=2,
    markerfacecolor="white",
    markersize=10,
)
plt.plot(
    [t_mem / t_flop, t_mem / t_flop],
    [1, 1e14],
    linewidth=2,
    linestyle=":",
    color="blue",
)


ax.set_xlabel("arithmetic intensity $q$")
ax.set_ylabel(r"performance $R [\operatorname{FLOPs}/\operatorname{s}]$")
plt.legend(loc="upper left", fontsize=14)
plt.savefig("roofline.svg", bbox_inches="tight")
