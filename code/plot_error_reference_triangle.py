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
plt.savefig("../lectures/figures/error_reference_triangle.png", bbox_inches="tight")
