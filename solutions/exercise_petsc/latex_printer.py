from petsc4py import PETSc
import numpy as np


def print_latex(x):
    if type(x) is PETSc.Mat:
        mat_dense = PETSc.Mat()
        x.convert("dense", mat_dense)
        a = mat_dense.getDenseArray()
        print(r"\begin{pmatrix}")
        for i in range(a.shape[0]):
            s = []
            for j in range(a.shape[1]):
                if abs(a[i, j]) < 1.0e-12:
                    s.append(r"\textcolor{lightgray}{0}")
                else:
                    s.append(f"{a[i,j]:4.1f}")
            print(" & ".join(s) + r"\\")
        print(r"\end{pmatrix}")

    elif type(x) is PETSc.Vec:
        a = x.getArray()
        print(r"\begin{pmatrix}")
        s = []
        for i in range(len(a)):
            if abs(a[i]) < 1.0e-12:
                s.append(r"\textcolor{lightgray}{0}")
            else:
                s.append(f"{a[i]:4.1f}")
        print(r" \\ ".join(s))
        print(r"\end{pmatrix}")
