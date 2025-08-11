"""Run numerical experiments (if necessary) and visualise results"""

import subprocess
import re
import json
import numpy as np
from matplotlib import pyplot as plt
import os


def run(n, ksp_type, pc_type):
    """Run numerical experiment for a specified setup

    Return tuple (niter, tolve) with number of iterations and
    solve time

    :arg n: problem size n
    :arg ksp_type: KSP type of solver
    :arg pc_type: PC type of solver
    """
    cmd = f"python linear_solve.py -n {n} -ksp_type {ksp_type} -pc_type {pc_type} -ksp_max_it 100000 -ksp_converged_reason"

    output = subprocess.run(cmd.split(), stdout=subprocess.PIPE).stdout.decode("utf-8")
    converged = True
    for line in output.split("\n"):
        m = re.match("number of iterations = *([0-9]+) *", line)
        if m:
            niter = int(m.group(1))
        m = re.match(r"time \[solve\] = *(.+) s", line)
        if m:
            tsolve = float(m.group(1))
        m = re.match(".*Linear solve did not converge.*", line)
        if m:
            converged = False
    if not converged:
        niter = -1
    return niter, tsolve


def plot_results(data, problemsizes, ylabel, filename, ylim, scale=1.0):
    """Plot the results for a specific data dictionary

    :arg data: dictionary such that data[ksp_type][pc_type] is a list
        of numbers with the quantity to plot
    :arg problemesizes: problem sizes
    :arg label: label to use on y-axis
    :arg filename: name of file to write to
    :arg ylim: tuple with range on y-axis
    :arg scale: scaling factor
    """
    plt.clf()
    fig, axs = plt.subplots(ncols=3, figsize=(12, 4))

    for j, ksp_type in enumerate(data.keys()):
        ax = axs[j]
        for pc_type in data[ksp_type].keys():
            Y = np.asarray(data[ksp_type][pc_type])
            ax.plot(
                problemsizes,
                scale * Y,
                linewidth=2,
                markersize=4,
                marker="o",
                label=f"{pc_type}",
            )
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.legend(loc="upper left")
        if j == 0:
            ax.set_ylabel(ylabel)
        ax.set_xlabel("problem size $n$")
        ax.set_title(f"{ksp_type}")
        ax.set_ylim(*ylim)
    plt.savefig(filename, bbox_inches="tight")


niter_dict = dict()
tsolve_dict = dict()
problemsizes = [16, 32, 64, 128, 256, 512, 1024, 2048]

# if data file does not exist generate it by running the code
if not os.path.exists("data.json"):
    for ksp_type in "richardson", "cg", "gmres":
        niter_dict[ksp_type] = {}
        tsolve_dict[ksp_type] = {}
        for pc_type in "jacobi", "ilu", "gamg":
            niter_dict[ksp_type][pc_type] = []
            tsolve_dict[ksp_type][pc_type] = []
            print(ksp_type, pc_type)
            for n in problemsizes:
                niter, tsolve = run(n, ksp_type, pc_type)
                print(f"{n:3d} : {niter:6d} {tsolve:12.8f} s")
                niter_dict[ksp_type][pc_type].append(niter)
                tsolve_dict[ksp_type][pc_type].append(tsolve)
            print()

    with open("data.json", "w") as f:
        json.dump(dict(niter=niter_dict, tsolve=tsolve_dict), f, indent=4)
else:
    with open("data.json", "r") as f:
        data = json.load(f)
        niter_dict = data["niter"]
        tsolve_dict = data["tsolve"]
        niter_dict = data["niter"]


plot_results(niter_dict, problemsizes, "number of iterations", "niter.svg", [1, 20000])
plot_results(
    tsolve_dict,
    problemsizes,
    r"solve time [$\mu \text{s}]",
    "tsolve.svg",
    [10, 300000],
    scale=1e6,
)
