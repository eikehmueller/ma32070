
<div align="center">
  <p style="font-size:32px;">MA32070 Exercise 5: Model solution</p>
</div>

----

*&#169; Eike Mueller, University of Bath 2025. These notes are copyright of Eike Mueller, University of Bath. They are provided exclusively for educational purposes at the University and are to be downloaded or copied for your private study only. Further distribution, e.g. by upload to external repositories, is prohibited. html generated with [pandoc](https://pandoc.org/) using [easy-pandoc-templates](https://github.com/ryangrose/easy-pandoc-templates) under the [GPL-3.0.1 license](https://github.com/ryangrose/easy-pandoc-templates?tab=GPL-3.0-1-ov-file#readme)*

----

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

The entire source code is available in [`linear_algebra.py`](linear_algebra.py).

# PETSc solver options

## Implementation
The source code is available in [`linear_solve.py`](linear_solve.py).

## Numerical experiments

### Gauss Seidel iteration
The problem $(D+L)\boldsymbol{z} = \boldsymbol{r}$ can be written in matrix form as

$$
\begin{pmatrix}
A_{00} & 0 & 0 & 0 &\dots & 0 \\
A_{10} & A_{11} & 0 &0 & & 0 \\
A_{20} & A_{21} & A_{22} &0 & & 0 \\
\vdots & & & & \ddots & \vdots\\ 
A_{n-1,0} & A_{n-1,1} & A_{n-1,2} & A_{n-1,3} & \dots & A_{n-1,n-1}
\end{pmatrix}
\begin{pmatrix}
z_0 \\ z_1 \\ z_2 \\\vdots \\ z_{n-1}
\end{pmatrix}
=
\begin{pmatrix}
r_0 \\ r_1 \\ r_2 \\\vdots \\ r_{n-1}
\end{pmatrix}
$$

The $i$-th equation is given by:

$$
\sum_{j=0}^{i-1} A_{ij}z_j + A_{ii} z_i = r_i
$$

Hence, if we know $z_0,z_1,\dots,z_{i-1}$ we can compute $z_i$ as

$$
z_i = \frac{1}{A_{ii}}\left(r_i-\sum_{j=0}^{i-1} A_{ij}z_j\right).
$$

Starting with $z_0 = r_0/A_{00}$, this can be used to compute $z_1,z_2,\dots,z_{n-1}$ recursively.

## Number of solver iterations for Jacobi and SOR
The following table shows the number of iterations required to solve the linear problem $A\boldsymbol{u}=\boldsymbol{b}$ to a tolerance of $10^{-6}$ for the Richardson iteration, preconditioned with Jacobi and SOR.

| problem size $n$ | number of iterations (Jacobi) | number of iterations (SOR) |
| ---:             | ---:                          | ---:                       |
| 16               |                         6 686 |                      1 538 |
| 32               |                        25 712 |                      6 147 |
| 64               |                      100 763  |                    26 117  |
| 128              |                       387 726 |                     89 239 |
| 256              |                     1 450 776 |                    387 421 |
| 512              |                     5 568 289 |                  1 478 886 |

For both preconditioners the number of iterations increases by a factor of approximately $4\times$ when the problem size is doubled. The absolute number of iterations is smaller for SOR than for Jacobi. However, as we will see below, a single application of the SOR preconditioner is significantly more expensive than applying the Jacobi method.

#### Number of solver iterations
![Number of solver iterations](niter.svg)
For the Richardson iteration, the number of solver iterations increases in proportion the the square of the problem size unless the GAMG preconditioner is used. As already observed above, the Jacobi preconditioner requires more iterations than SOR. For the ICC method the number of iterations also increases with problem size but not quite as strongly as for Jacobi and SOR. The ILU preconditioner does not work in this setup.

The CG solver reduces the number of iterations dramatically, but it still grows for Jacobi and SOR. For ICC and GAMG the number of iterations is independent of the problem size. Again ILU does not work as a preconditioner and this is not very surprising since CG requires a symmetric preconditioner, which ILU is not.

ICC and GAMG perform very similarly when used with GMRES instead of CG. For this solver, SOR and ILU can be used, but only for smaller problem sizes. In both cases, the number of iterations increases with problem size. Using the Jacobi method results in a rapid increase in the number of iterations.

As far as the number of iterations in concerned, the best setup appears to be a combination of CG or GMRES with either ICC or GAMG.
#### Total solution time
![Total solution time](tsolve.svg)
in all cases the total solution time increases with problem size, and this can be attributed to a combination of the increase in number of iterations (previous section) and the growing cost per iteration (next section). The use of CG or GMRES dramatically improves overall performance, in particular when combined with the ICC preconditioner. Although both ICC and GMRES require comparable numbers of iterations, ICC is faster overall since an application of a single GAMG solve is relatively expensive. Observe also that when used with the CG solver, the Jacobi iteration can be competitive, in particular for smaller problems size. The reason for this is that the Jacobi preconditioner only requires division by the diagonal of $A$, which is a very cheap operation.

#### Time per solver iteration
![Time per solver iteration](titer.svg)
For sufficiently large $n$, the cost per iteration grows in direct proportion to the problem size. While for a given preconditioner the time per iteration depends on the iterative solver, this dependence is relatively weak. This implies that the majority of time is spent in the preconditioner application. As expected, a single application of simple preconditioners such as Jacobi and SOR is relatively cheap, but this is counteracted by the large number of iterations required to achieve convergence. The more powerful preconditioners (GAMG and ICC) are significantly more expensive but result in the smallest overall solution time since they require only a small number of iterations to converge.