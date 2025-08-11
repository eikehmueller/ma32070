import sys
import petsc4py

petsc4py.init(sys.argv)

from petsc4py import PETSc
import numpy as np

n = 512
h_sq = 1 / n**2

col_indices = list(
    np.asarray(
        [[(row + col) % n for col in range(-1, 2)] for row in range(n)]
    ).flatten()
)
row_start = list(np.arange(0, n * 3, 3)) + [3 * n]
values = n * [-1, 2 + 1 * h_sq, -1]

A = PETSc.Mat().createAIJWithArrays(size=(n, n), csr=(row_start, col_indices, values))

ksp = PETSc.KSP().create()
ksp.setOperators(A)
ksp.setFromOptions()
rng = np.random.default_rng(seed=1241773)
b = PETSc.Vec().createWithArray(rng.normal(size=n))
u = PETSc.Vec().createWithArray(np.zeros(n))
ksp.solve(b, u)

A.convert("dense")
mat = A.getDenseArray()
nrow, ncol = mat.shape

epsilon = 1e-12
for i in range(nrow):

    s = " & ".join(
        [
            f"{x:3.0f}" if abs(x) > epsilon else r"\textcolor{lightgray}{0}"
            for x in A[i, :]
        ]
    )
    s += r"\\"
    # print(s)
