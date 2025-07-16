# Sparse matrix representations
The stiffness matrices we obtain from out finite element discretisation contain a lot of zero entries. Consider, for example, the $81\times 81$ matrix that is obtained for a piecewise linear discretisation on a $8\times 8$ grid:

![Stiffness matrix](figures/stiffness_matrix.png)

Of the $81\times 81 = 6561$ entries of this matrix, only $n_{\text{nz}}=497$ or $7.6\%$ are nonzero, which corresponds to an average number of around $\overline{n}_{\text{nz}} = 6.14$ nonzeros per row. For large $n$, the number of non-zeros per row will remain roughly the same, leading to an even poorer fill-ratio.

Clearly, it is very inefficient to store all these zero entries if only $\mathcal{O}(n)$ entries are in fact required to encode the data stored in the matrix.

## Compressed Sparse Row storage
A matrix $A$ with $n_{\text{nz}}\ll n$ nonzero entries is often called *sparse*. To store sparse matrices, we can proceed as follows:

1. Store all non-zero entries in a long array $V$ of length $n_{\text{nz}}$, going throw the matrix row by row
2. For each non-zero entry also store the column index in an array $J$ of the same length
 
We could now also store the corresponding row-indices in an array $I$ of the same length. With this, it would then be possible to reconstruct all non-zero entries of $A$:

#### Algorithm: Reconstruction of matrix 
1. Set $A\gets 0$
2. for $\ell=0,1,2,\dots,n_{\text{nz}}$ **do**
3. $~~~~$ Set $A_{I_\ell,J_\ell} \gets V_{\ell}$
4. **end do**

However, there is a more efficient way of doing this: Since the arrays $V$ and $J$ are constructed by going through the matrix row by row, we only need to keep track of the positions where a new row starts. This can be encoded as follows:

3. Store an array $R$ of length $n+1$ such that $R_i$ describes the index in $V$, $J$ where a new row starts. For convenience, we also store $R_{n} = n_{\text{nz}}$.

The resulting storage format, consisting of the arrays $V$ (values), $J$ (column indices) and $R$ (row pointers) is known as Compressed Sparse Row storage (CSR)

#### Algorithm: Reconstruction of matrix in CSR
1. Set $A\gets 0$
2. Set $\ell\gets 0$
3. for $i=0,1,2,\dots,n-1$ **do**
4. $~~~~$ for $j=R_i,R_i+1,\dots,R_{i+1}-1$ **do**
5. $~~~~~~~~$ Set $A_{i,J_\ell} \gets V_{\ell}$
6. $~~~~~~~~$ Increment $\ell\gets \ell+1$
7. $~~~~$ **end do**
8. **end do**

#### Example
Consider the following $5\times 5$ matrix with $n_{\text{nz}}=11$ non-zero entries:

$$
\begin{pmatrix}
1 & 2 & 0 & 8 & 0 \\
4 & 6 & 0 & 0 & 0 \\
0 & 2 & 8 & 0 & 9 \\
0 & 0 & 0 & 0 & 0 \\
0 & 3 & 1 & 0 & 7
\end{pmatrix}
$$

We have the following arrays:

1. Values: $V=[1,2,8,4,6,2,8,9,3,1,7]$
2. Column indices: $J=[0,1,3,0,1,1,2,4,1,2,4]$
3. Row pointers: $R=[0,3,6,8,8,11]$

Note that one of the rows contains only zero entries.

### Matrix-vector multiplication

#### Algorithm: Matrix-vector multiplication $y = y + Ax$ in CSR storage
1. Set $\ell\gets 0$
2. for $i=0,1,2,\dots,n-1$ **do**
3. $~~~~$ for $j=R_i,R_i+1,\dots,R_{i+1}-1$ **do**
4. $~~~~~~~~$ Set $y_i \gets y_i + V_{\ell} x_{J_{\ell}}$
5. $~~~~~~~~$ Increment $\ell\gets \ell+1$
6. $~~~~$ **end do**
7. **end do**

## PETSc implementation
To implement matrices in the CSR storage format, we use the [PETSc](https://petsc.org) library, more explicitly, the [petsc4py](https://petsc.org/release/petsc4py/) Python interface. After installation, this can be imported as follows:
```python
from petsc4py import PETSc
```
We can now create an (empty) matrix with

```python
matrix = PETSc.Mat()
```

To create the $5\times 5$ matrix above we first need to set up the sparsity structure, i.e. the arrays $J$ (`col_indices`) and $R$ (`row_start`):
```python
n_row = 5
n_col = 5

col_indices = [0, 1, 3, 0, 1, 1, 2, 4, 1, 2, 4]
row_start = [0, 3, 5, 8, 8, 11]

matrix.createAIJ((n_row, n_col), csr=(row_start, col_indices))
```
We can now insert values, for example we might want to set $A_{0,3} = 8$, as highlighted in red here:
$$
\begin{pmatrix}
1 & 2 & 0 & \textcolor{red}{8} & 0 \\
4 & 6 & 0 & 0 & 0 \\
0 & 2 & 8 & 0 & 9 \\
0 & 0 & 0 & 0 & 0 \\
0 & 3 & 1 & 0 & 7
\end{pmatrix}
$$
This can be done by calling the `setValue()` method:
```python
# Set 
row = 0
col = 3
value = 8
matrix.setValue(row, col, value)
```
Note that this method has an optional parameter `addv`, for `addv=True` the value will be added to an already existing value.

We can also set blocks of value. For example, we might want to set the $2\times 2$ block in the upper left corner, as highlighhted in red here:
$$
\begin{pmatrix}
\textcolor{red}{1} & \textcolor{red}{2} & 0 & 8 & 0 \\
\textcolor{red}{4} & \textcolor{red}{6} & 0 & 0 & 0 \\
0 & 2 & 8 & 0 & 9 \\
0 & 0 & 0 & 0 & 0 \\
0 & 3 & 1 & 0 & 7
\end{pmatrix}
$$
For this, we need to specify the rows and columns in the target matrix as follows:
```python
rows = [0, 1]
cols = [0, 1]
local_matrix = np.asarray([1, 2, 4, 6])
matrix.setValues(rows, cols, local_matrix)
```
Blocks do not have to be contiguous. We could, for example, set the 6 non-zero values highlighted in red here:
$$
\begin{pmatrix}
1 & 2 & 0 & 8 & 0 \\
4 & 6 & 0 & 0 & 0 \\
0 & \textcolor{red}{2} & \textcolor{red}{8} & 0 & \textcolor{red}{9} \\
0 & 0 & 0 & 0 & 0 \\
0 & \textcolor{red}{3} & \textcolor{red}{1} & 0 & \textcolor{red}{7}
\end{pmatrix}
$$
The indices of these values are described by the tensor product $(2,4)\times(1,2,4)$, and hence we need to do this:
```python
rows = [2, 4]
cols = [1, 2, 4]
local_matrix = np.asarray([2, 8, 9, 3, 1, 7])
matrix.setValues(rows, cols, local_matrix)
```
Finally, before we can use the matrix, we need to assemble it:
```python
matrix.assemble()
```
For more information also see the [`petsc4py.Mat` documentation](https://petsc.org/release/petsc4py/reference/petsc4py.PETSc.Mat.html).

### Matrix-vector multiplication
PETSc also provides a vector class. For example, to create the vector
$$
v = \begin{pmatrix}
8\\0\\9\\4\\5
\end{pmatrix}
$$
we can do this:
```python
v = PETSc.Vec()
v.createWithArray([8, 0, 9, 4, 5])
```
We can now multiply the matrix that we created above with this vector to compute $w=Av$:
```python
w = matrix @ v
```
To print the vector we need to first extract the underlying array with the `getArray()` method:
```python
print(w.getArray())
```

