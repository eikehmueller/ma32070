import numpy as np

epsilon = 4.0e-16

theta = np.pi / 8
D = np.asarray([[1 + epsilon, 1], [1, 1]], dtype=np.float64)
R = np.asarray([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])

A = np.asarray(
    [
        [
            1 + np.sqrt(2) / 2 + (2 + np.sqrt(2)) / 4 * epsilon,
            (1 - epsilon / 2) * np.sqrt(2) / 2,
        ],
        [
            (1 - epsilon / 2) * np.sqrt(2) / 2,
            1 - np.sqrt(2) / 2 + (2 - np.sqrt(2)) / 4 * epsilon,
        ],
    ],
    # dtype=np.float32,
)
A_star = R @ D @ R.T
print(A - A_star)
# A = A_star
# b = R @ np.asarray([1, 1])
b = np.asarray(
    [
        (np.sqrt(2 + np.sqrt(2)) + np.sqrt(2 - np.sqrt(2))) / 2,
        (np.sqrt(2 + np.sqrt(2)) - np.sqrt(2 - np.sqrt(2))) / 2,
    ]
)

x = np.linalg.solve(A, b)
x_true = np.asarray([np.sqrt(2 - np.sqrt(2)) / 2, np.sqrt(2 + np.sqrt(2)) / 2])


print(f"b      = {b[0]:18.16f} {b[1]:18.16f}")
print(f"x      = {x[0]:18.16f} {x[1]:18.16f}")
print(f"x_true = {x_true[0]:18.16f} {x_true[1]:18.16f}")
error_nrm = np.linalg.norm(x - x_true) / np.linalg.norm(x_true)
print(f"error norm = {error_nrm:8.4e}")
print()
C = (2 - np.sqrt(2)) / 4
print(f"{C:18.16}")
# print(np.linalg.eigvals(A))
