import numpy as np

rng = np.random.default_rng(seed=361847)

u = rng.normal(size=4)
v = rng.normal(size=4)

A = rng.normal(size=(3, 4))
B = rng.normal(size=(4, 3))

T = rng.normal(size=(4, 3, 2))
S = rng.normal(size=(2, 3, 5))
Q = rng.normal(size=(5, 3, 2, 2))

# A u
print("A u = ")
print(A @ u)
print(np.dot(A, u))
print(np.einsum("ij,j->i", A, u))
print()

# A B
print("A B")
print(A @ B)
print(np.dot(A, B))
print(np.einsum("ik,kj->ij", A, B))
print()

# u x v
print("u x v")
print(np.einsum("i,j->ij", u, v))
print(np.tensordot(u, v, axes=0))
print()

# trace(A B^T)
print("trace(A B^T)")
print(np.einsum("ij,ji->", A, B))
print(np.sum(A * B.T))
print()

print("T S Q")
print(np.einsum("aji,bjk,kjii->ab", T, S, Q))
print()
