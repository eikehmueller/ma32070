<div align="center">
  <p style="font-size:32px;">MA32070 Exercise 1: Model solution</p>
</div>

----

*&#169; Eike Mueller, University of Bath 2025. These notes are copyright of Eike Mueller, University of Bath. They are provided exclusively for educational purposes at the University and are to be downloaded or copied for your private study only. Further distribution, e.g. by upload to external repositories, is prohibited. html generated with [pandoc](https://pandoc.org/) using [easy-pandoc-templates](https://github.com/ryangrose/easy-pandoc-templates) under the [GPL-3.0.1 license](https://github.com/ryangrose/easy-pandoc-templates?tab=GPL-3.0-1-ov-file#readme)*

----

# Implementation

The class `CubicElement` is a subclass of the abstract `FiniteElement` class, which needs to be imported from [fem/finiteelement.py]()https://github.com/eikehmueller/finiteelements/blob/main/src/fem/finiteelement.py:

```Python
from fem.finiteelement import FiniteElement

class CubicElement(FiniteElement):  
```
## Class methods
We need to make sure that concrete implementations are provided for all abstract methods in `FiniteElement`. For this, we first implement a method which computes the Vandermonde matrix for a given set of points.

### Vandermonde matrix
For the given $\nu=10$ monomials $\{\theta_\ell(x)\}_{\ell=0}^{9} = \{1,x_0,x_1,x_0^2,x_0x_1,x_1^2,x_0^3,x_0^2x_1,x_0x_1^2,x_1^3\}$ and $n$ points $\zeta^{(0)},\zeta^{(1)},\dots,\zeta^{(n-1)}$ the Vandermonde matrix has shape $(n,10)$ and the following entries:

$$
V = V(\{\zeta^{(r)}\}_{r=0}^{n-1}) = \begin{pmatrix}
  1 & \zeta^{(0)}_0 & \zeta^{(0)}_1 & (\zeta^{(0)}_0)^2 & \zeta^{(0)}_0 \zeta^{(0)}_1 & (\zeta^{(0)}_1)^2 & (\zeta^{(0)}_0)^3 & (\zeta^{(0)}_0)^2 \zeta^{(0)}_1 & \zeta^{(0)}_0 (\zeta^{(0)}_1)^2 & (\zeta^{(0)}_1)^3\\
  1 & \zeta^{(1)}_0 & \zeta^{(1)}_1 & (\zeta^{(1)}_0)^2 & \zeta^{(1)}_0 \zeta^{(1)}_1 & (\zeta^{(1)}_1)^2 & (\zeta^{(1)}_0)^3 & (\zeta^{(1)}_0)^2 \zeta^{(1)}_1 & \zeta^{(1)}_0 (\zeta^{(1)}_1)^2 & (\zeta^{(1)}_1)^3\\
  \vdots & & & & & \vdots & & & & \vdots\\
  1 & \zeta^{(n-1)}_0 & \zeta^{(n-1)}_1 & (\zeta^{(n-1)}_0)^2 & \zeta^{(n-1)}_0 \zeta^{(n-1)}_1 & (\zeta^{(n-1)}_1)^2 & (\zeta^{(n-1)}_0)^3 & (\zeta^{(n-1)}_0)^2 \zeta^{(n-1)}_1 & \zeta^{(n-1)}_0 (\zeta^{(n-1)}_1)^2 & (\zeta^{(n-1)}_1)^3\\
\end{pmatrix}\qquad(\ddagger)
$$

Hence, for `grad=False` we return the matrix `mat` which can be constructed like so:

```Python
mat = np.empty([npoints, 10])
mat[:, 0] = 1
mat[:, 1] = zeta[..., 0]
mat[:, 2] = zeta[..., 1]
mat[:, 3] = zeta[..., 0] ** 2
mat[:, 4] = zeta[..., 0] * zeta[..., 1]
mat[:, 5] = zeta[..., 1] ** 2
mat[:, 6] = zeta[..., 0] ** 3
mat[:, 7] = zeta[..., 0] ** 2 * zeta[..., 1]
mat[:, 8] = zeta[..., 0] * zeta[..., 1] ** 2
mat[:, 9] = zeta[..., 1] ** 3
```

For `grad=True` the method should return the rank-3 tensor $V^\partial$ of shape $(n,10,2)$ with

$$
V^\partial_{rma} = \frac{\partial \theta_m}{\partial x_a}(\zeta^{(r)})
$$

We can interpret $V^\partial$ as two $n\times 10$ matrices $V^{\partial_0}=\frac{\partial V}{\partial x_0}$, $V^{\partial_1}=\frac{\partial V}{\partial x_1}$ obtained by taking suitable pertial derivatives of the Vandermonde matrix $V$ defined in $(\dagger)$ above. Explicitly these matrices are given by

$$
V^{\partial_0}  = \begin{pmatrix}
  0 &  1 & 0 & 2\zeta^{(0)}_0 & \zeta^{(0)}_1  & 0 &3(\zeta^{(0)}_0)^2 & 2\zeta^{(0)}_0 \zeta^{(0)}_1 & (\zeta^{(0)}_1)^2 & 0\\
  0 &  1 & 0 & 2\zeta^{(1)}_0 & \zeta^{(1)}_1  & 0 &3(\zeta^{(1)}_0)^2 & 2\zeta^{(1)}_0 \zeta^{(1)}_1 & (\zeta^{(1)}_1)^2 & 0\\
  \vdots & & & & & \vdots & & & & \vdots\\
  0 &  1 & 0 & 2\zeta^{(n-1)}_0 & \zeta^{(n-1)}_1  & 0 &3(\zeta^{(n-1)}_0)^2 & 2\zeta^{(n-1)}_0 \zeta^{(n-1)}_1 & (\zeta^{(n-1)}_1)^2 & 0\\
\end{pmatrix}\\[4ex]
V^{\partial_1} = \begin{pmatrix}
  0 & 0 & 1 & 0 & \zeta^{(0)}_0 & 2\zeta^{(0)}_1 & 0 & (\zeta^{(0)}_0)^2 & 2\zeta^{(0)}_0\zeta^{(0)}_1 & 3 (\zeta^{(0)}_1)^2\\
  0 & 0 & 1 & 0 & \zeta^{(1)}_0 & 2\zeta^{(1)}_1 & 0 & (\zeta^{(1)}_0)^2 & 2\zeta^{(1)}_0\zeta^{(1)}_1 & 3 (\zeta^{(1)}_1)^2\\
  \vdots & & & & & \vdots & & & & \vdots\\
  0 & 0 & 1 & 0 & \zeta^{(n-1)}_0 & 2\zeta^{(n-1)}_1 & 0 & (\zeta^{(n-1)}_0)^2 & 2\zeta^{(n-1)}_0\zeta^{(n-1)}_1 & 3 (\zeta^{(n-1)}_1)^2
\end{pmatrix}
$$
In the code, the matrix $V^{\partial}$ is constructed with
```Python
mat = np.zeros([npoints, 10, 2])
mat[:, 1, 0] = 1
mat[:, 3, 0] = 2 * zeta[..., 0]
mat[:, 4, 0] = zeta[..., 1]
mat[:, 6, 0] = 3 * zeta[..., 0] ** 2
mat[:, 7, 0] = 2 * zeta[..., 0] * zeta[..., 1]
mat[:, 8, 0] = zeta[..., 1] ** 2

mat[:, 2, 1] = 1
mat[:, 4, 1] = zeta[..., 0]
mat[:, 5, 1] = 2 * zeta[..., 1]
mat[:, 7, 1] = zeta[..., 0] ** 2
mat[:, 8, 1] = 2 * zeta[..., 0] * zeta[..., 1]
mat[:, 9, 1] = 3 * zeta[..., 1] ** 2
```

### Tabulation of basis functions
To tabulate the basis functions at some points $\{\zeta^{(r)}\}_{r=0}^{n-1}$, we need to evaluate the Vandermonde matrix at these points and multiply it by the coefficient matrix $C$ to obtain $T=V(\{\zeta^{(r)}\}_{r=0}^{n-1})C$. This is readily done with the `@` operator for matrix-matrix multiplication:

```Python
_zeta = np.asarray(zeta)
mat = self._vandermonde_matrix(_zeta,grad=False) @ self._coefficients
```

This assumes that `zeta` is passed as an array of shape $(n,2)$. If `zeta` is a single point, i.e. a two-dimensional vector, we first use [`np.expand_dims()`](https://numpy.org/doc/stable/reference/generated/numpy.expand_dims.html) to convert it to an array of shape $(1,2)$ which can be passed to the `_vandermonde_matrix()` method. The resulting matrix `mat` will have shape $(1,nu)$ and we use [`np.squeeze()`](https://numpy.org/doc/stable/reference/generated/numpy.squeeze.html) to turn it into a $\nu$-dimensional vector. Instead of using an if-statement to distinguish the two cases, we can write this more compactly as

```Python
mat = np.squeeze(
    self._vandermonde_matrix(
        np.expand_dims(_zeta, axis=list(range(2 - _zeta.ndim))),
        grad=False,
    )
    @ self._coefficients
)
```

### Tabulation of gradients of basis functions
To compute $T^\partial$ with $T^{\partial}_{r\ell a} = \sum_{m=0}^{\nu-1} V^\partial_{rma}(\{\zeta^{(r)}\}_{r=0}^{n-1}) C_{m\ell}$ we use [np.einsum()](https://numpy.org/doc/stable/reference/generated/numpy.einsum.html). Again, `np.expand_dims()` and `np.squeeze()` are used to deal with the special case where `zeta` is a vector.

```Python
mat = np.squeeze(
    np.einsum(
        "imk,mj->ijk",
        self._vandermonde_matrix(
            np.expand_dims(_zeta, axis=list(range(2 - _zeta.ndim))),
            grad=True,
        ),
        self._coefficients,
    )
)
```

### Tabulation of dofs
Tabulation of the dofs is achieved by evaluating the passed function $\widehat{f}$ at the nodal points $\xi$:

```Python
def tabulate_dofs(self, fhat):
    return fhat(self._nodal_points)
```
The nodal points $\{\xi^{(\ell)}\}_{\ell=0}^{\nu-1}$ can be initialised when the object is created by setting 

```Python
self._nodal_points = np.asarray(
            [
                [0, 0],
                [1, 0],
                [0, 1],
                [2 / 3, 1 / 3],
                [1 / 3, 2 / 3],
                [0, 2 / 3],
                [0, 1 / 3],
                [1 / 3, 0],
                [2 / 3, 0],
                [1 / 3, 1 / 3],
            ]
        )
```
in the `__init__()` method. In the same method we also construct the coefficient matrix $C = \left(V(\{\xi^{(\ell)}\}_{\ell=0}^{\nu-1})\right)^{-1}$ (stored in `self._oefficient`), which is given by the inverse of the Vandermonde matrix $V(\{\xi^{(\ell)}\}_{\ell=0}^{\nu-1})$ evaluated at the nodal points :

```Python
vandermonde_matrix = self._vandermonde_matrix(self._nodal_points)
self._coefficients = np.linalg.inv(vandermonde_matrix)
```

### Class properties
The implementation of the abstract properties `ndof_per_interior`, `ndof_per_facet` and `ndof_per_vertex`, which should return 1, 2 and 1 respectively, is straightforward.


## Complete code
The full implementation of the `CubicElement` class can be found in [cubicelement.py](cubicelement.py). The corresponding tests are contained in [test_cubicelement.py](test_cubicelement.py).