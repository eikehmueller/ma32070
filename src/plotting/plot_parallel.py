import numpy as np
from matplotlib import pyplot as plt
import glob
import re


def t_matmat(n, p):
    """Theoretical time for square matrix-matrix multiplication

    :arg n: matrix size
    :arg p: number of processors
    """
    t_flop = 1.77e-11
    t_mem = 2.74e-10
    t_lat = 1.34e-6
    t_word = 2.7e-9
    return (2 * n**3 * t_flop + 3 * n**2 * t_mem) / p + n**2 * t_word + p * t_lat


def extract_measured_data(directory, problemsizes=None, overlap=False):
    filenames = glob.glob(
        f"{directory}/runtime_" + ("overlap_" if overlap else "") + "p*.txt"
    )
    pss = set()
    processors = set()
    for filename in filenames:
        m = re.match(r".*_p([0-9]+)_n([0-9]+)\.txt", filename)
        if m:
            processors.add(int(m.group(1)))
            pss.add(int(m.group(2)))
    pss = problemsizes if problemsizes else sorted(list(pss))
    processors = np.asarray(sorted(list(processors)))
    T_dict = {}
    for n in pss:
        data = []
        for p in processors:
            filename = (
                f"{directory}/runtime_"
                + ("overlap_" if overlap else "")
                + f"p{p}_n{n}.txt"
            )
            with open(filename, "r", encoding="utf8") as f:
                for line in f.readlines():
                    m = re.match(" *elapsed time *= *(.+) *s", line)
                    if m:
                        data.append(float(m.group(1)))
        T_dict[n] = np.asarray(data)
    return processors, T_dict


def plot_performance(P, T_dict, filename):
    problemsizes = T_dict.keys()
    plt.clf()
    fig, axs = plt.subplots(1, 3, figsize=(4, 3))
    fig.subplots_adjust(left=0.2, right=3.4)
    marker = ["o", "s", "v", "^", "<", ">", "D"]
    for j, n in enumerate(problemsizes):
        T = T_dict[n]
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
            T_dict[n][0] / T,
            linewidth=2,
            marker=marker[j],
            markersize=6,
            markeredgewidth=2,
            markerfacecolor="white",
            label=f"$n={n}$",
        )
        axs[2].plot(
            P,
            T_dict[n][0] / (T * P),
            linewidth=2,
            marker=marker[j],
            markersize=6,
            markeredgewidth=2,
            markerfacecolor="white",
        )
    axs[1].plot(P, P, linewidth=2, color="black", linestyle="--")
    for ax in axs:
        ax.set_xlim(0.5, 2**10)
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
    plt.savefig(filename, bbox_inches="tight")


P = 2 ** np.arange(10)
problemsizes = [256, 1024, 4096, 16384]
T_dict = {n: t_matmat(n, P) for n in problemsizes}
plot_performance(P, T_dict, "parallel_scaling_theory.pdf")
P, T_dict = extract_measured_data(
    "output_parallel", problemsizes=problemsizes, overlap=False
)
plot_performance(P, T_dict, "parallel_scaling_measured.pdf")
P, T_dict = extract_measured_data(
    "output_parallel", problemsizes=problemsizes, overlap=True
)
plot_performance(P, T_dict, "parallel_scaling_measured_overlap.pdf")
