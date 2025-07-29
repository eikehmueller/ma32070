import numpy as np


def to_latex(matrix, color_row=np.inf):
    _matrix = np.array(matrix)
    if _matrix.ndim == 1:
        _matrix = np.expand_dims(_matrix, -1)
    nrow, ncol = _matrix.shape
    s = ""
    s += r"\begin{pmatrix}" + "\n"
    for i in range(nrow):
        row = [
            f"{x:.4g}" if abs(x) > 1e-12 else "\\textcolor{lightgray}{0}"
            for x in _matrix[i, :]
        ]
        if i >= color_row:
            for j in range(color_row, ncol):
                if i == j and i == color_row:
                    color = "red"
                elif i == color_row:
                    color = "green"
                else:
                    color = "blue"
                row[j] = f"\\textcolor{{{color}}}{{{row[j]}}}"
        s += " & ".join(row)
        if i < nrow - 1:
            s += r"\\"
        s += "\n"
    s += r"\end{pmatrix}"
    return s


n = 5
rng = np.random.default_rng(seed=312847)

A = (2 / n * rng.uniform(size=(n, n), low=0, high=1) + np.eye(n)).round(decimals=4)
b = rng.uniform(size=n, low=0, high=1).round(decimals=2)

A_in = np.matrix(A)
b_in = np.matrix(b)

u = np.linalg.solve(A, b)

print("% Initial")
print("$$")
print("\\begin{aligned}")
print("    A &= ", to_latex(A, color_row=0), "\\\\")
print("    b &= ", to_latex(b, color_row=0))
print("\\end{aligned}")
print("$$")


for k in range(0, n - 1):
    for i in range(k + 1, n):
        pivot = A[i, k] / A[k, k]
        A[i, k:] -= pivot * A[k, k:]
        b[i] -= pivot * b[k]
    print(f"% After stage {k+1}:")
    print("$$")
    print("\\begin{aligned}")
    print("    A &= ", to_latex(A, color_row=k + 1), "\\\\")
    print("    b &= ", to_latex(b, color_row=k + 1))
    print("\\end{aligned}")
    print("$$")
    print()

print("solution:")
print("    u = ", to_latex(u))

assert np.allclose(A @ u, b, rtol=1e-12)
assert np.allclose(A_in @ u, b_in, rtol=1e-12)
