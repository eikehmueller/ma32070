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
        1: 7.5080e-01,
        2: 2.0408e-03,
        3: 2.1487e-04,
        4: 3.9261e-05,
        5: 1.0675e-06,
        6: 3.6658e-07,
        7: 3.2341e-09,
        8: 1.9690e-09,
        9: 8.5817e-12,
        10: 7.7725e-12,
        11: 3.3446e-14,
        12: 2.5368e-14,
        13: 6.6124e-14,
        14: 8.6731e-12,
        15: 3.1712e-10,
        16: 2.7028e-08,
        17: 1.2541e-06,
        18: 6.0339e-05,
        19: 5.3898e-02,
        20: 8.1112e-01,
    },
    "single": {
        1: 7.5080e-01,
        2: 2.0407e-03,
        3: 2.1489e-04,
        4: 3.9258e-05,
        5: 1.0625e-06,
        6: 3.6783e-07,
        7: 3.1128e-09,
        8: 2.6805e-09,
        9: 4.0985e-10,
        10: 5.4061e-07,
        11: 4.0023e-06,
        12: 3.2557e-05,
        13: 3.9902e-04,
        14: 1.2133e-03,
        15: 3.9817e-02,
        16: 2.2864e-03,
        17: 1.3432e00,
        18: 3.3853e-01,
        19: 1.0613e00,
        20: 9.6616e-01,
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
ax.set_ylabel(r"error $\|u-u_{\text{exact}}\|^2_{L_2}/\|u_{\text{excat}}\|^2_{L_2}$")
ax.set_xticks(range(1, 21), labels=[str(j) for j in range(1, 21)], rotation=15)
plt.legend(loc="lower left")
plt.savefig("error_reference_triangle.png", bbox_inches="tight")
