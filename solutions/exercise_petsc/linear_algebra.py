from petsc4py import PETSc


def print_as_dense(x):
    """Print object

    If x is a sparse matrix, first convert it to a dense matrix.

    :arg x: object to print, either PETSc.Mat or PETSc.Vec
    """
    if type(x) is PETSc.Mat:
        mat_dense = PETSc.Mat()
        x.convert("dense", mat_dense)
        print(mat_dense.getDenseArray())
    elif type(x) is PETSc.Vec:
        print(x.getArray())


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
print("matrix A")
print_as_dense(A)
print()

# Create 4x4 matrix B
B = PETSc.Mat()
col_indices = [0, 1, 0, 2, 3]
row_start = [0, 2, 3, 4, 5]
values = [-0.7, 2.5, 8.7, 3.2, 12.0]
B.createAIJWithArrays((n_row, n_col), (row_start, col_indices, values))
B.assemble()
print("matrix B")
print_as_dense(B)
print()

# Create the vectors v and w
v = PETSc.Vec()
v.createWithArray([7.3, -0.7, 0, 3.2])
print("vector v")
print_as_dense(v)
print()

w = PETSc.Vec()
w.createWithArray([0, 0, 0.3, 2.8])
print("vector w")
print_as_dense(w)
print()

# Compute matrix-matrix product A B
C = A.matMult(B)
print("matrix A B (using matMult())")
print_as_dense(C)
print()
C = A @ B
print("matrix A B (using @ operator)")
print_as_dense(C)
print()

# Compute matrix-matrix product A B^T
C = A.matTransposeMult(B)
print("matrix A B^T (using matTransposeMult() )")
print_as_dense(C)
print()
C = A @ B.transpose()
print("matrix A B^T (using @ operator )")
print_as_dense(C)
print()

# Compute matrix-matrix product A^T B
C = A.transposeMatMult(B)
print("matrix A B^T (using transposeMatMult() )")
print_as_dense(C)
print()
C = A @ B.transpose()
print("matrix A B^T (using @ operator)")
print_as_dense(C)
print()

# Compute sum A + B^T
C = A + B.transpose()
print("matrix A + B^T")
print_as_dense(C)
print()

# Compute matrix-vector product Av + w
r = PETSc.Vec()
r.createSeq(n_row)
A.multAdd(v, w, r)
print("vector A v + w (using multAdd() )")
print_as_dense(r)
print()
r = A @ v + w
print("vector A v + w (using @ operator)")
print_as_dense(r)
print()

# Compute matrix-vector product B^T w
r = PETSc.Vec()
r.createSeq(n_row)
B.multTranspose(w, r)
print("vector B^T w (using multTranspose() )")
print_as_dense(r)
print()
r = B.transpose() @ w
print("vector B^T w (using @ operator)")
print_as_dense(r)
print()
