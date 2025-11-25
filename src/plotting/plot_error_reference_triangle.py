"""# width of peak
sigma = 0.5
# location of peak
x0 = np.asarray([0.6, 0.25])
# Coefficient of diffusion term
kappa = 0.9
# Coefficient of zero order term
omega = 0.4
"""

import numpy as np
from matplotlib import pyplot as plt

rel_error_norm_squared = {
    "double": {
        1: 8.2010e-03,
        2: 7.5380e-04,
        3: 7.0303e-05,
        4: 7.9261e-06,
        5: 1.9660e-07,
        6: 3.6353e-08,
        7: 3.1679e-10,
        8: 9.6407e-11,
        9: 3.7619e-13,
        10: 1.6767e-13,
        11: 5.7372e-16,
        12: 3.2589e-16,
        13: 4.0343e-15,
        14: 2.8875e-13,
        15: 9.5789e-11,
        16: 9.0484e-09,
        17: 1.9771e-07,
        18: 1.6864e-06,
        19: 8.3504e-03,
        20: 2.4018e-01,
    },
    "single": {
        1: 8.2010e-03,
        2: 7.5379e-04,
        3: 7.0303e-05,
        4: 7.9260e-06,
        5: 1.9660e-07,
        6: 3.6820e-08,
        7: 4.8390e-10,
        8: 2.8462e-10,
        9: 3.4271e-11,
        10: 9.9478e-08,
        11: 9.5775e-07,
        12: 1.2327e-05,
        13: 1.1919e-04,
        14: 2.2606e-04,
        15: 2.9481e-03,
        16: 5.1202e00,
        17: 5.3173e-01,
        18: 1.9845e00,
        19: 1.5743e-01,
        20: 2.4328e-01,
    },
}

condition_number = {
    1: 20.50000000000002,
    2: 132.29717461240412,
    3: 490.9053014732179,
    4: 1353.5364303993297,
    5: 3540.0042929515253,
    6: 9620.383349672371,
    7: 27688.708734093445,
    8: 84905.26933463394,
    9: 274533.2695542549,
    10: 926297.2982223562,
    11: 3221906.0005568536,
    12: 11462411.711292647,
    13: 41442441.35018524,
    14: 151715174.92217693,
    15: 560713932.5821387,
    16: 2088876768.1209428,
    17: 7830494218.344533,
    18: 29218311962.73611,
    19: 79471534108.5321,
    20: 118793192661.2284,
}

plt.clf()
color = {"double": "blue", "single": "red"}
for precision, data in rel_error_norm_squared.items():
    X = np.asarray(list(data.keys()))
    Y = np.asarray(list(data.values()))
    plt.plot(
        X,
        Y,
        linewidth=2,
        color=color[precision],
        marker="o",
        markersize=4,
        label=precision,
    )
ax = plt.gca()
ax.set_yscale("log")
ax.set_xlabel("polynomial degree $p$")
ax.set_ylabel(r"error $\|u-u_{\text{exact}}\|^2_{L_2(\widehat{K})}$")
ax.set_xticks(range(1, 21), labels=[str(j) for j in range(1, 21)], rotation=15)
plt.legend(loc="lower left")
plt.savefig("../../lectures/figures/error_reference_triangle.png", bbox_inches="tight")

plt.clf()
color = {"double": "blue", "single": "red"}
X = np.asarray(list(condition_number.keys()))
Y = np.asarray(list(condition_number.values()))
epsilon = {"single": 1.19e-7, "double": 2.22e-16}
N_dof = (X + 1) * (X + 2) / 2
plt.plot(
    X,
    N_dof * Y,
    linewidth=2,
    color="blue",
    marker="o",
    markersize=4,
)

ax = plt.gca()
ax.set_yscale("log")
ax.set_xlabel("polynomial degree $p$")
ax.set_ylabel(r"$n_{\text{dof}}\cdot \text{cond}(A^{(h)})$")
ax.set_xticks(range(1, 21), labels=[str(j) for j in range(1, 21)], rotation=15)
plt.savefig(
    "../../lectures/figures/conditioning_reference_triangle.png", bbox_inches="tight"
)
