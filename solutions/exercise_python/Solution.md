# Solution

* [numerical_integration.py](numerical_integration.py)
* [driver.py](driver.py)
* [linear_algebra_numpy.py](linear_algebra_numpy.py)

# Linear algebra with PETSc
The matrices can be created either with `createAIJ()` followed by calls to `setValue()` (in this case, `assemble()` needs to be called to complete the assembly process ) or in one go with `createAIJWithArrays()`. Below, we use the two different approaches for $A$ and $B$ respectively.

```Python
# Number of rows and columns
n_row = 4
n_col = 4

# Create 4x4 matrix A
A = PETSc.Mat()
col_indices = [0, 1, 1, 3, 2, 0, 3]
row_start = [0, 2, 4, 5, 7]
A.createAIJ((n_row, n_col), csr=(row_start, col_indices))
A.setValue(0, 0, 1.7)
A.setValue(0, 1, 2.3)
A.setValue(1, 1, -3.4)
A.setValue(1, 3, 4.5)
A.setValue(2, 2, 8.6)
A.setValue(3, 0, -1.3)
A.setValue(3, 3, 1.2)
A.assemble()

# Create 4x4 matrix B
B = PETSc.Mat()
col_indices = [0, 1, 0, 2, 3]
row_start = [0, 2, 3, 4, 5]
values = [-0.7, 2.5, 8.7, 3.2, 12.0]
B.createAIJWithArrays((n_row, n_col), (row_start, col_indices, values))
```
The vectors can be created with `createWithArray()`, for example

```Python
v = PETSc.Vec()
v.createWithArray([7.3, -0.7, 0, 3.2])
```

The results of the linear algebra operations are as follows:

The product $AB$ can be computed either with `C = A.matMult(B)` or by using the `@` operator as `C = A @ B`.
$$
A B =\begin{pmatrix}
18.8 &  4.2 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0}\\
-29.6 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 54.0\\
\textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 27.5 & \textcolor{lightgray}{0}\\
 0.9 & -3.2 & \textcolor{lightgray}{0} & 14.4\\
\end{pmatrix}
$$

The product $AB^\top$ can be computed with `C = A.matTransposeMult(B)` or as `C = A @ B.transpose()`.
$$ 
A B^T = \begin{pmatrix}
 4.6 & 14.8 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0}\\
-8.5 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 54.0\\
\textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 27.5 & \textcolor{lightgray}{0}\\
 0.9 & -11.3 & \textcolor{lightgray}{0} & 14.4\\
\end{pmatrix}
$$
Similarly, the product $A^\top B$ can be computed with `C = A.transposeMatMult(B)` or as `C = A.transpose() @ B`.
$$
A^T B = \begin{pmatrix}
-1.2 &  4.2 & \textcolor{lightgray}{0} & -15.6\\
-31.2 &  5.8 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0}\\
\textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 27.5 & \textcolor{lightgray}{0}\\
39.1 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 14.4\\
\end{pmatrix}
$$
To add the transpose of matrix the matrix $B$ to the matrix $A$ we can use `A + B.transpose()`.
$$
A + B^T = \begin{pmatrix}
 1.0 & 11.0 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0}\\
 2.5 & -3.4 & \textcolor{lightgray}{0} &  4.5\\
\textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 11.8 & \textcolor{lightgray}{0}\\
-1.3 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 13.2\\
\end{pmatrix}
$$
The matrix-vector product $A\boldsymbol{v} + \boldsymbol{w}$ can be computed by constructing a vector `r` with and then calling the `multAdd()` method of `A`:
```
r = PETSc.Vec()
r.createSeq(n_row)
A.multAdd(v, w, r)
```
Alternatively, we can also simply use `r = A @ v + w`, in this case it is not necessary to create `r` as a PETSc vector first.
$$
A \boldsymbol{v} + \boldsymbol{w} = \begin{pmatrix}
10.8 \\ 16.8 \\  0.3 \\ -2.9
\end{pmatrix}
$$
Similarly, to compute $B^T\boldsymbol{w}$ we can use either `B.multTranspose(w, r)` or simply `B.transpose() @ w`.
$$
B^T \boldsymbol{w} = \begin{pmatrix}
\textcolor{lightgray}{0} \\ \textcolor{lightgray}{0} \\  1.0 \\ 33.6
\end{pmatrix}
$$

The entire source code is available in [`linear_algebra_petsc.py`](linear_algebra_petsc.py).