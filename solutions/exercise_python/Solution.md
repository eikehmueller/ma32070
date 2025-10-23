# Numerical integration

## Implementation of `evaluate()`
The `evaluate()` method of the `NumericalIntegration` base class can be implemented with a for loop:

```Python
def evaluate(self, f):
    """Numerically approximate the integral int_a^b f(x) dx

    For this, loop over the subintervals [x_j,x_{j+1}], approximate
    the integral in each subinterval and sum the result

    :arg f: function to integrate
    """
    h = (self._b - self._a) / self._n
    integral = 0.0
    for j in range(self._n):
        integral += self._integrate(f, self._a + j * h, self._a + (j + 1) * h)
    return integral
```

Alternatively, we could use a list comprehension to construct a list `Is` which contains the integrals over the sub-intervals $[I(x_0,x_1), I(x_1,x_2), \dots, I(x_{n-1},x_n)]$ and then sum all entries of the list with `np.sum()` as follows:

```Python
Is = [ self._integrate(f, self._a + j * h, self._a + (j + 1) * h) for j in range(self._n) ]
integral = float(np.sum(Is))
```

We convert the result to the native Python `float` type. In fact, this approach is not ideal, since for large $n$ the list `Is` is constructed and stored explicitly in memory, while we only need to inspect one element at a time and add it to the sum. To avoid this, one could use a [generator expression](https://docs.python.org/3/reference/expressions.html#generator-expressions). This is achieved by replacing the square brackets `[`,`]` in the definition of `Is` with round brackets `(`, `)`:

```Python
Is = ( self._integrate(f, self._a + j * h, self._a + (j + 1) * h) for j in range(self._n) )
```
In this case, `Is` is a "recipe" for constructing the values of $I(x_j,x_{j+1})$ which can be consumed by `np.sum()`. See [the documentation](https://wiki.python.org/moin/Generators) for more on Python generators.

## Implementation of midpoint- and Simpsons-rule
The classes `MidpointRule` and `SimpsonsRule` are both derived from `NumericalIntegration`. Both require the implementation of the `_integrate()` method. 

```Python
class MidpointRule(NumericalIntegration):
    """Numerical integration with the midpoint rule"""
    def __init__(self, interval, n):
        """Initialise instance

        :arg interval: interval [a,b]
        :arg n: number of subintervals
        """
        super().__init__(interval, n, order=2)

    def _integrate(self, f, x_m, x_p):
        """Approximate int_{x_-}^{x_+} f(x)dx by midpoint rule

        (x_+ - x_-) * f((x_+ + x_-)/2)

        :arg f: function to integrate
        :arg x_m: lower bound x_-
        :arg x_p: upper bound x_+
        """
        return (x_p - x_m) * f((x_m + x_p) / 2)

class SimpsonsRule(NumericalIntegration):
    """Numerical integration with Simpson's rule"""
    def __init__(self, interval, n):
        """Initialise instance

        :arg interval: interval [a,b]
        :arg n: number of subintervals
        """
        super().__init__(interval, n, order=4)

    def _integrate(self, f, x_m, x_p):
        """Approximate int_{x_-}^{x_+} f(x)dx by Simpson's rule

        (x_+ - x_-)/6 * ( f(x_-) + 4*f((x_+ + x_-)/2) + f(x_+) )

        :arg f: function to integrate
        :arg x_m: lower bound x_-
        :arg x_p: upper bound x_+
        """
        return (x_p - x_m) / 6 * (f(x_m) + 4 * f((x_m + x_p) / 2) + f(x_p))
```

The `order` property is set by passing it to the initialiser of the base class when calling `super().__init__()`.

The complete source code can be found in [numerical_integration.py](numerical_integration.py).

## Implementation of main program
The `results` dictionary can be construced by appending to the lists for each integrator in a for-loop like this:
```Python
results = {"MidpointRule": [], "SimpsonsRule": []}
for n in [4, 8, 16, 32]:
    integrator_midpoint = MidpointRule([0, 1], n)
    results["MidpointRule"].append(
        integrator_midpoint.evaluate(functools.partial(f, alpha=alpha))
    )
for n in [4, 8, 16, 32]:
    integrator_simpson = SimpsonsRule([0, 1], n)
    results["SimpsonsRule"].append(
        integrator_simpson.evaluate(functools.partial(f, alpha=alpha))
    )
```

However, it is more elegant to use list- and [dictionary](https://docs.python.org/3/tutorial/datastructures.html#dictionaries) comprehensions like this:
```Python
results = {
    Integrator.__name__: [
        Integrator([0, 1], n).evaluate(functools.partial(f, alpha=alpha))
        for n in [4, 8, 16, 32]
    ]
    for Integrator in [MidpointRule, SimpsonsRule]
}
```

#### Print statement
To understand the code for printing out the results, namely
```Python
for integrator, integrals in results.items():
    print(
        f"{integrator}: "
        + ", ".join([f"{abs(x-exact_result):8.4e}" for x in integrals])
    )
```
observe that the outer for-loop iterates over the keys (the `integrator` strings) and values (the lists `integrals`) of the `results` dictionary.

In the print-statement
```Python
f"{integrator}: "
```
constructs a string with the name of the integrator followed by `: `. The code 
```Python
[f"{abs(I-I_exact):8.4e}" for I in integrals]
```
uses a list comprehension to construct a list of strings, each of which is obtained by computing $\left|I^{\text{(integrator)}}(n) - I\right|$ and formatting the result as a floating point number of 8 charactes in total and 4 places after the decimal point. The resulting list of strings is passed to the [str.join()](https://docs.python.org/3/library/stdtypes.html#str.join) method which concatenates the strings and separates them by `, `:
```Python
", ".join([f"{abs(I-I_exact):8.4e}" for I in integrals])
```

The complete source code can be found in [driver.py](driver.py).

### Order of convergence
To empirically compute the order of convergence, observe that for sufficiently small $h$ we have that

$$
I(n) \approx I + C h^\mu \qquad\text{with $h=(b-a)/n$.}
$$

It is easy to see that

$$
\mu \approx \log_2\left( \frac{I(n)-I(2n)}{I(2n)-I(4n)}\right)
$$

Crucially, the exact value (which might be unknown) and the constants $C$ cancel out in the ratio. The following code can be used to compute and print $\mu$ for different values of $n$:

```Python
for integrator, integrals in results.items():
    print(f"{integrator}: ")
    mu = [
        np.log2(
            (integrals[j] - integrals[j + 1]) / (integrals[j + 1] - integrals[j + 2])
        )
        for j in range(len(integrals) - 2)
    ]
    print(
        "  mu = ",
        ", ".join([f"{z:8.4f}" for z in mu]),
    )
```
We obtain the following results, which confirm that the integrators are indeed of order $2$ and $4$ respectively:

| $n$   | Midpoint rule | Simpson's rule |
| :---: | :---:         | :---:          |
| $4$   | $1.9996$      | $3.9997$       |   
| $8$   | $1.9999$      | $3.9999$       |

# Linear algebra with numpy
The matrix-vector product $A\boldsymbol{u}$ can be computed in (at least) three different ways:

1. `A @ u`
2. `np.dot(A, u)`
3. `np.einsum("ij,j->i", A, u)`

The matrix-matrix product $AB$ can be computed with

1. `A @ B`
2. `np.dot(A, B)`
3. `np.einsum("ik,kj->ij", A, B)`

The outer product $\boldsymbol{u}\otimes \boldsymbol{v}$ can be computed with

1. `np.einsum("i,j->ij", u, v)`
2. `np.tensordot(u, v, axes=0)`

The trace $\text{trace}(AB)$ can be computed with

1. `np.sum(A * B.T)`
2. `np.einsum("ij,ji->", A, B)`
3. `np.linalg.trace(A@B)`

Finally, the specified triple-product of $T$, $S$ and $Q$ can be computed with
`np.einsum("aji,bjk,kjii->ab", T, S, Q)`.

The complete source code can be found in [linear_algebra_numpy.py](linear_algebra_numpy.py).

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

The product $AB$ can be computed either with [`matMult()`](https://petsc.org/release/petsc4py/reference/petsc4py.PETSc.Mat.html#petsc4py.PETSc.Mat.matMult) as `C = A.matMult(B)` or by using the `@` operator as `C = A @ B`.
$$
A B =\begin{pmatrix}
18.8 &  4.2 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0}\\
-29.6 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 54.0\\
\textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 27.5 & \textcolor{lightgray}{0}\\
 0.9 & -3.2 & \textcolor{lightgray}{0} & 14.4\\
\end{pmatrix}
$$

The product $AB^\top$ can be computed with [`matTransposeMult()`](https://petsc.org/release/petsc4py/reference/petsc4py.PETSc.Mat.html#petsc4py.PETSc.Mat.matTransposeMult) as `C = A.matTransposeMult(B)` or by using the `@` operator and [`createTranspose()`](https://petsc.org/release/petsc4py/reference/petsc4py.PETSc.Mat.html#petsc4py.PETSc.Mat.createTranspose) as `C = A @ PETSc.Mat().createTranspose(B)`.

$$ 
A B^T = \begin{pmatrix}
 4.6 & 14.8 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0}\\
-8.5 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 54.0\\
\textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 27.5 & \textcolor{lightgray}{0}\\
 0.9 & -11.3 & \textcolor{lightgray}{0} & 14.4\\
\end{pmatrix}
$$

Observe that `PETSc.Mat().createTranspose(B)` creates a new matrix to represent $B^\top$, whereas [`transpose()`](https://petsc.org/release/petsc4py/reference/petsc4py.PETSc.Mat.html#petsc4py.PETSc.Mat.transpose) performs the transition in-place, i.e. overwrites $B$ with its transpose. Hence, `C = A @ B.transpose()` would also work, but after executing this line the variable `B` will store $B^\top$ so care needs to be taken if this matrix is used for further calculations.

Similarly, the product $A^\top B$ can be computed with [`transposeMatMult()`](https://petsc.org/release/petsc4py/reference/petsc4py.PETSc.Mat.html#petsc4py.PETSc.Mat.transposeMatMult) as `C = A.transposeMatMult(B)` or as `C = PETSc.Mat().createTranspose(A) @ B`.
$$
A^T B = \begin{pmatrix}
-1.2 &  4.2 & \textcolor{lightgray}{0} & -15.6\\
-31.2 &  5.8 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0}\\
\textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 27.5 & \textcolor{lightgray}{0}\\
39.1 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 14.4\\
\end{pmatrix}
$$
To add the transpose of matrix the matrix $B$ to the matrix $A$ we can use `A + PETSc.Mat().createTranspose(B)`.
$$
A + B^T = \begin{pmatrix}
 1.0 & 11.0 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0}\\
 2.5 & -3.4 & \textcolor{lightgray}{0} &  4.5\\
\textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 11.8 & \textcolor{lightgray}{0}\\
-1.3 & \textcolor{lightgray}{0} & \textcolor{lightgray}{0} & 13.2\\
\end{pmatrix}
$$
The matrix-vector product $A\boldsymbol{v} + \boldsymbol{w}$ can be computed by constructing a vector `r` with [`createSeq()`](https://petsc.org/release/petsc4py/reference/petsc4py.PETSc.Vec.html#petsc4py.PETSc.Vec.createSeq) and then calling the [`multAdd()`](https://petsc.org/release/petsc4py/reference/petsc4py.PETSc.Mat.html#petsc4py.PETSc.Mat.multAdd) method of `A`:
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
Similarly, to compute $B^T\boldsymbol{w}$ we can use either [`multTranspose()`](https://petsc.org/release/petsc4py/reference/petsc4py.PETSc.Mat.html#petsc4py.PETSc.Mat.multTranspose) as `B.multTranspose(w, r)` or simply `PETSc.Mat().createTranspose(B) @ w`.
$$
B^T \boldsymbol{w} = \begin{pmatrix}
\textcolor{lightgray}{0} \\ \textcolor{lightgray}{0} \\  1.0 \\ 33.6
\end{pmatrix}
$$

The entire source code is available in [linear_algebra_petsc.py](linear_algebra_petsc.py).