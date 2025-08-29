"""Plot performance of parallel matrix-matrix multiplication"""

import numpy as np
from matplotlib import pyplot as plt
import glob
import re


def t_matmat(n, p):
    """Theoretical time for square matrix-matrix multiplication

    :arg n: matrix size
    :arg p: number of processors
    """
    t_flop = 0.2444e-10
    t_mem = 0.9241e-10
    t_lat = 1.1896e-06
    t_word = 9.9629e-10
    return (2 * n**3 * t_flop + 3 * n**2 * t_mem) / p + n**2 * t_word + p * t_lat


def extract_measured_data(directory, overlap=False):
    """Extract measured runtimes from files of the form

    runtime_p{P}_n{N}.txt or runtime_overlap_p{P}_n{N}.txt

    where P is the number of processors and N the problem size. A typical
    file runtime_p8_n256.txt will have the following content:

        Running on 8 processors.
        n =      256
        m =      256
        r =      256
        overlap = False
        nocomm  = False
        error norm = 5.7071e-16
        elapsed time = 5.01e-04 s

    :arg directory: directory which contains the output files
    :arg overlap: process results for runs with overlap of comm and comms
    """
    filenames = glob.glob(
        f"{directory}/runtime_" + ("overlap_" if overlap else "") + "p*.txt"
    )
    problemsizes = set()
    processors = set()
    for filename in filenames:
        m = re.match(r".*_p([0-9]+)_n([0-9]+)\.txt", filename)
        if m:
            processors.add(int(m.group(1)))
            problemsizes.add(int(m.group(2)))
    problemsizes = sorted(list(problemsizes))
    processors = np.asarray(sorted(list(processors)))
    results = {}
    for n in problemsizes:
        runtime = []
        nproc = []
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
                        nproc.append(p)
                        runtime.append(float(m.group(1)))
        results[n] = [np.asarray(nproc), np.asarray(runtime)]
    return results


def plot_performance(P, results, filename):
    plt.clf()
    fig, axs = plt.subplots(1, 3, figsize=(4, 3))
    fig.subplots_adjust(left=0.2, right=3.4)
    marker = ["o", "s", "v", "^", "<", ">", "D"]
    for j, (n, data) in enumerate(results.items()):
        P, T = data
        p_time = axs[0].plot(
            P,
            T,
            linewidth=2,
            marker=marker[j],
            markersize=6,
            markeredgewidth=2,
            markerfacecolor="white",
        )
        axs[0].plot(
            P,
            T[0] / P,
            linewidth=2,
            markersize=0,
            linestyle="--",
            color=p_time[0].get_color(),
        )
        axs[1].plot(
            P,
            T[0] / T,
            linewidth=2,
            marker=marker[j],
            markersize=6,
            markeredgewidth=2,
            markerfacecolor="white",
            label=f"$n={n}$",
        )
        axs[2].plot(
            P,
            T[0] / (T * P),
            linewidth=2,
            marker=marker[j],
            markersize=6,
            markeredgewidth=2,
            markerfacecolor="white",
        )
    axs[1].plot([1, 2**10], [1, 2**10], linewidth=2, color="black", linestyle="--")
    for ax in axs:
        ax.set_xlim(0.5, 2**7)
        ax.set_xscale("log")
        ax.set_xlabel("number of processors")
    axs[0].set_yscale("log")
    axs[0].set_ylabel("runtime $T_p$ [s]")
    axs[0].set_ylim(1e-5, 400)
    axs[1].set_ylabel("speedup $T_1/T_p$")
    axs[1].set_yscale("log")
    axs[1].set_ylim(0.5, 2**7)
    axs[1].legend(loc="upper left")
    axs[2].set_ylabel(r"parallel efficiency $T_1/(p\cdot T_p)$")
    axs[2].set_yticks([0, 0.25, 0.5, 0.75, 1])
    axs[2].set_ylim(0, 1.1)
    axs[2].set_yticklabels([r"0%", r"25%", r"50%", r"75%", r"100%"])
    plt.savefig(filename, bbox_inches="tight")


def fit_pingpong(data_filename, filename):
    n = []
    t = []
    with open(data_filename, "r", encoding="utf8") as f:
        for line in f.readlines():
            m = re.match(
                r"^ *Python message +size *= *([0-9]+) *t *= *([0-9\.]+) *micro +seconds *$",
                line,
            )
            if m:
                n.append(int(m.group(1)))
                t.append(1e-6 * float(m.group(2)))
    n = np.asarray(n)
    t = np.asarray(t)
    plt.clf()
    p_fit = np.polyfit(n[:], t[:], deg=1)
    N = np.arange(n[0], n[-1], 1)
    plt.plot(n, t, linewidth=2, markersize=6, marker="o")
    plt.plot(N, p_fit[0] * N + p_fit[1], linewidth=2, color="black", linestyle="--")
    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.savefig(filename, bbox_inches="tight")
    t_lat = p_fit[1]
    t_mem = p_fit[0]
    print(f"t_lat = {t_lat:8.4e} s")
    print(f"t_mem = {t_mem:8.4e} s")


extension = "svg"

fit_pingpong("output_parallel/performance.txt", f"pingpong.{extension}")
P = 2 ** np.arange(7)
problemsizes = [256, 1024, 4096, 16384]
results_theory = {n: [P, t_matmat(n, P)] for n in problemsizes}
plot_performance(P, results_theory, f"parallel_scaling_theory.{extension}")
results_measured = extract_measured_data("output_parallel", overlap=False)
results_measured = {n: results_measured[n] for n in problemsizes}
plot_performance(P, results_measured, f"parallel_scaling_measured.{extension}")
results_measured_overlap = extract_measured_data("output_parallel", overlap=True)
results_measured_overlap = {n: results_measured_overlap[n] for n in problemsizes}
plot_performance(
    P, results_measured_overlap, f"parallel_scaling_measured_overlap.{extension}"
)
