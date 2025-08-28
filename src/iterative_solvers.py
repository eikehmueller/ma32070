import sys
import petsc4py

petsc4py.init(sys.argv)

from petsc4py import PETSc

row_start = [0, 3, 5, 6, 8, 9]
col_indices = [0, 1, 3, 0, 1, 2, 0, 3, 4]
values = [10.2, 0.8, -2.1, 0.8, 6.7, 6.4, -2.1, 7.2, 9.8]

A = PETSc.Mat().createAIJWithArrays(
    (5, 5),
    (
        row_start,
        col_indices,
        values,
    ),
)
A.assemble()
A_dense = PETSc.Mat()
A_dense = A.convert("dense")
A_numpy = A_dense.getDenseArray()
import numpy as np

print(A_numpy)
# assert False

b = PETSc.Vec().createWithArray([8.1, 0, 9.3, -4.3, 5.2])
u = PETSc.Vec().createSeq(b.size)

ksp = PETSc.KSP().create()
ksp.setOperators(A)
ksp.setFromOptions()
ksp.solve(b, u)
