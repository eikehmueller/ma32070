# Finite Elements
We start by solving the weak form of the PDE for a special case, namely a domain consisting of a single triangle. By doing this, we will develop the fundamental concepts and techniques of finite element analysis and discuss their implementation in Python. As we will see later, the solution of the PDE on more complicated domains can be reduced to this case.
## Triangular reference element
Let us consider a very simple domain $K=\Omega$ which consists of the unit triangle with vertices $v_0=(0,0)$, $v_1=(1,0)$ and $v_2=(0,1)$. We label the edges (or facets) in a counter-clockwise fashion as $F_0 = \overrightarrow{v_1v_2}$, $F_1 = \overrightarrow{v_2v_0}$ and $F_2 = \overrightarrow{v_0v_1}$:

![reference triangle](figures/reference_triangle.svg)

In the following we will also refer to this as the *reference triangle*.

Recall that the finite element approach starts with the choice of a suitable function space $V$. For this, consider the space of bi-variate polynomials of degree $p$ on $K$:

$$
\mathcal{P}_p(K) = \{q:q(x) = \sum_{\substack{\alpha_0,\alpha_1\\\alpha_0+\alpha_1\le p}} a_{\alpha_0,\alpha_1} x_0^{\alpha_0}x_1^{\alpha_1},\;a_{\alpha_0,\alpha_1}\in\mathbb{R}\}\subset H^1(K)
$$

The space $\mathcal{P}_p(K)$  is spanned by $d_p = {p+2 \choose 2} = \frac{1}{2}(p+2)(p+1)$ basis functions $\{\phi_j(x)\}_{j=0}^{d_p-1}$. These can be chosen to be the monomials $\{1,x_0,x_1,x_0^2,x_0x_1,x_1^2,\dots$\}, but a better choice is to pick [Lagrange polynomials](https://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html). This will later allow us to construct $H^1(\Omega)$ functions on a mesh that consists of little triangles by "glueing together" the functions on neighbouring triangles. To construct Lagrange polynomials, we choose $d_p$ points $\{\xi^{(j)}\}_{j=0}^{d_p-1}$ in $K$ and define $\phi_j(x)\in\mathcal{P}_p(K)$ such that
$$
\phi_j(\xi^{(k)}) = \delta_{jk} = \begin{cases}
    1 & \text{for $j=k$}\\
    0 & \text{otherwise}.
\end{cases}
$$
A possible choice of points is given by
$$
\{\xi^{(j)}\}_{j=0}^{d_p-1} = \left\{\left(\frac{j_0}{p},\frac{j_1}{p}\right) \quad \text{for $0\le j_0\le j_1 \le p$}\right\}.
$$

Note that this choice of points is not optimal for higher polynomial degrees $p$ in the sense that it can lead to numerical instabilities. For all problems we consider in this course this does not matter, however. 

We order these points (and the associated basis functions $\phi_j(x)$) as follows:

* points associated with the three vertices $v_0$, $v_1$, $v_2$ (in this order); there is $\nu_{\text{vertex}}=1$ point per vertex
* points associated with the facets $F_0$, $F_1$, $F_2$ (in this order); there are $\nu_{\text{facet}}=p-1$ points per facet and on each facet these points are ordered according to the arrows in the figure above
* points associated with the interior $K^0$; there are $\nu_{\text{interior}}=\frac{1}{2}(p-1)(p-2)$ points of this type

This is illustrated in the following figure, which also shows the ordering of the nodes for $p=1,2,3,4$:

![Nodes of Lagrane elements](figures/lagrange_nodes.svg)

The associated finite elements are also known as **Lagrange elements**.

### Examples
#### Linear finite element
For $p=1$ we obtain the (bi-)linear finite element with the following basis three functions:
$$
\begin{aligned}
\phi_0(x) &= 1-x_0-y_0\\
\phi_1(x) &= x_0\\
\phi_2(x) &= x_1
\end{aligned}
$$
#### Quadratic finite element
For $=2$ there are six basis functions, three associated with vertices
$$
\begin{aligned}
\phi_0(x,y) &= (1-x-y)(1-2x-2y),\\
\phi_1(x,y) &= x(2x-1),\\
\phi_2(x,y) &= y(2y-1),
\end{aligned}
$$
and three associated with facets
$$
\begin{aligned}
\phi_3(x,y) &= 4xy,\\
\phi_4(x,y) &= 4y(1-x-y),\\
\phi_5(x,y) &= 4x(1-x-y).
\end{aligned}
$$
These functions are visualised in the following figure (red arrows indicate gradients):

![Basis functions for quadratic finite element](figures/quadratic_element.png)


## Formal definition of finite elements
It turns out that it is advantageous to define finite elements in a more general sense. Mirroring this more abstract mathematical definition in the Python code will help us to structure the code in a sensible way that will allow its easy adaptation to specific cases. For this we first need to introduce the notion of the dual $V^*$ for a function space $V$.
### Dual spaces
Consider a domain $K$ and the space $\mathcal{V}=\mathcal{V}(K)$ of real-valued functions $w:K\rightarrow \mathbb{R}$ on $K$. A linear *functional* $\ell$ maps a function $w\in \mathcal{V}$ to a real value such that

$$
\ell(c_1 w^{(1)}+c_2 w^{(2)}) = c_1\ell(w^{(1)})+c_2 \ell(w^{(2)}) \qquad\text{for all $c_1,c_2\in\mathbb{R}$, $w^{(1)}, w^{(2)} \in \mathcal{V}$}
$$

The space of all linear functionals on $\mathcal{V}$ is called the **dual space** $\mathcal{V}^*$.

#### Examples
Let $K\subset \mathbb{R}^2$ and $\mathcal{V}=H^1(K)$ be the space of functions with a square integrable first derivative. Then the following $\ell$ are linear functionals:
* point evaluation: $\ell(w) := w(\xi)$ for some point $\xi\in K$
* differentiation: $\ell(w) := \frac{\partial w}{\partial x_0}$
* integration: $\ell(w) := \int_K f(x)w(x)$ for some function $f(x)\in L_2(K)$

### Ciarlet's definition of the finite element
This now leads to the following definition (see [Logg2012]): a finite element is a triple $(K,\mathcal{V},\mathcal{L})$ which consists of
* the **domain** $K$
* the **function space** $\mathcal{V}=\mathcal{V}(K)$ of real-valued functions on $K$,
* the **degrees of freedom** (or **nodes**) $\mathcal{L} = \{\ell_j\}_{j=0}^{d-1}$, which is a basis for $\mathcal{V}^*$, the dual of $\mathcal{V}$

Crucially, we define the finite element by choosing a basis of the *dual* space $\mathcal{V}^*$. However, we can always construct a so-called *nodal* basis $\{\phi_j\}_{j=0}^{d-1}$ of $\mathcal{V}$ by requiring that
$$
\ell_j (\phi_k) = \delta_{jk} \qquad\text{for all $j,k=0,1,\dots,d-1$}.
$$
In the following we will assume that $K$ is the reference triangle introduced above, unless specified otherwise.

#### Examples
The **polynomial Lagrange element** we described above is a special case of this with 
* $\mathcal{V} = \mathcal{P}_p(K)$, the space of bi-variate polynomials of degree $p$
* $\ell_j: \ell_j(w) = w(\xi^{(j)})$ the point evaluation at the nodal points $\xi^{(j)}$

An alternative choice for the nodes would have been to define for some point $\widehat{\xi}\in K$:
$$
\ell_j (w) = \frac{\partial^{j_a}w}{\partial x_0^{j_b} \partial x_1^{j_a-j_b}}(\widehat{\xi}) \qquad\text{for $0\le j_b \le j_a\le p$ and $j=\frac{1}{2}j_a(j_a-1) + j_b$}
$$

The **Argyris finite element** (see Section 3.7.1 in [Logg2012]) is given by
* $\mathcal{V} = \mathcal{P}_5(K)$, the space of quintic bi-variate polynomials
* the 21 nodes defined as follows:
  - $\ell_i(w) = w(v_j)$ (evaluation at each vertex $v_i$ $\Rightarrow$ 3 nodes)
  - $\ell_{3+2i+j}(w) = \frac{\partial w}{\partial x_j}(v_i)$ (two gradient evaluations at each vertex $\Rightarrow$ 6 nodes)
  - $\ell_{9+3i+2j+k}(w) = \frac{\partial^2 w}{\partial x_j \partial x_k}(v_i)$ with $0\le j\le k\le 1$ (Hessian evaluation at each vertex $\Rightarrow$ 9 nodes)
  - $\ell_{18+j}(w) = n_i\cdot \nabla w(m_i)$ (normal derivative evaluation at the midpoints $m_i$ of each facet $F_i$ $\Rightarrow$ 3 nodes)

Note that the Argyris element and the quintic Lagrange element only differ in the choice of nodes. It turns out that the Argyris allows the construction of function spaces that have a bounded second derivative.

### Node numbering
As we will see later, it is crucial to establish a consistent ordering of the degrees of freedom. For this, assume that each node is associated with a topological entity of the reference triangle $K$. These entities are
* the vertices $v_0$, $v_1$, $v_2$ in this order
* the facets $F_0$, $F_1$, $F_2$ in this order
* the interior $K^0$ of $K$
We further assume that $0\le \nu_{\text{vertex}}$ nodes are associated with each vertex, $0\le \nu_{\text{facet}}$ nodes are associated with each facet and $0\le \nu_{\text{interior}}$ nodes are associated with the interior $K^0$. Then obviously
$$
d = 3( \nu_{\text{vertex}}+\nu_{\text{facet}})+\nu_{\text{interior}}.
$$
Let $\ell_k^{(E_i)}$ be the $k$-th node associated with topological entity $E\in \{v_0,v_1,v_2,F_0,F_1,F_2,K^0\}$. Then we arrange the unknowns $\{\ell_0,\dots,\ell_{d-1}\}$ in the following order:

$\{\ell_0^{(v_0)},\dots,\ell_{\nu_{\text{vertex}}-1}^{(v_0)},
\ell_0^{(v_1)},\dots,\ell_{\nu_{\text{vertex}}-1}^{(v_1)},
\ell_0^{(v_2)},\dots,\ell_{\nu_{\text{vertex}}-1}^{(v_2)},
\ell_0^{(F_0)},\dots,\ell_{\nu_{\text{facet}}-1}^{(F_0)},
\ell_0^{(F_1)},\dots,\ell_{\nu_{\text{facet}}-1}^{(F_1)},
\ell_0^{(F_2)},\dots,\ell_{\nu_{\text{facet}}-1}^{(F_2)},
\ell_0^{(K^0)},\dots,\ell_{\nu_{\text{interior}}-1}^{(K^0)}
\}$

In other words $\ell_{j=\mu_{\text{dof}}(E,i,k)} = \ell_k^{(E_i)}$ with the indirection map
$$
\mu_{\text{dof}}(E,i,k) = \begin{cases}
i\cdot \nu_{\text{vertex}} + k & \text{if $E=v_i$}\\
3\nu_{\text{vertex}} + i\cdot \nu_{\text{facet}} + k & \text{if $E=F_i$}\\
3(\nu_{\text{vertex}} + \nu_{\text{facet}}) + i & \text{if $E=K^0$}
\end{cases}
$$
This is illustrated for the polynomial Lagrange element in the figure above.

## Vandermonde matrix
Having picked the nodes, how can we construct the nodal basis functions $\{\phi_k(x)\}_{k=0}^{d-1}$ for a given set of nodes $\{\ell_j\}_{j=0}^{d-1}$? For this, assume that we know some set of basis functions $\{b_i(x)\}_{i=0}^{d-1}$ of $\mathcal{V}$. For the Lagrange elements, these could for example be the monomials $1,x_0,x_1,x_0^2,x_0x_1,x_1^2,\dots$. Hence, for each $k=0,1,\dots,d-1$ we can write

$$
\phi_k(x) = \sum_{i=0}^{d-1} c_i^{(k)} b_i(x)
$$

for some coefficients $c_i^{(k)}$. Since per definition $\{\phi_k\}_{k=0}^{d-1}$ is a nodal basis and $\ell_j$ are linear functionals we know that

$$
\delta_{jk} = \ell_j(\phi_k) = \sum_{i=0}^{d-1} \underbrace{c_i^{(k)}}_{C_{ik}} \underbrace{\ell_j(b_i)}_{V_{ji}}.
$$

If we define the $d\times d$ matrices $V$, $C$ with $V_{ji} := \ell_j(b_i)$ and $C_{ik}:=c_i^{(k)}$, then this equation can be written in matrix form as
$$
VC = \mathbb{I}\quad \Leftrightarrow \quad C = V^{-1}
$$
with $\mathbb{I}$ the $d\times d$ identity matrix. In other words, we can obtain the coefficients $c_i^{(k)}$ by inverting the matrix $V$. For the Lagrange element, where $\ell_j(w) = w(\xi^{(j)})$ are nodal evaluations, the matrix $V$ is the Vandermonde matrix:
$$
V = V(\{\xi^{(j)}\}_{j=0}^{d_p-1}) = \begin{pmatrix}
1 & \xi^{(0)}_0 & \xi^{(0)}_1  & (\xi^{(0)}_0)^2 & \xi^{(0)}_0 \xi^{(0)}_1 & (\xi^{(0)}_1)^2 & \dots \\[1ex]
1 & \xi^{(1)}_0 & \xi^{(1)}_1 & (\xi^{(1)}_0)^2 & \xi^{(1)}_0 \xi^{(1)}_1 & (\xi^{(1)}_1)^2 & \dots \\[1ex]
1 & \xi^{(2)}_0 & \xi^{(2)}_1 & (\xi^{(2)}_0)^2 & \xi^{(2)}_0 \xi^{(2)}_1 & (\xi^{(2)}_1)^2 & \dots \\[1ex]
\vdots & \vdots & \vdots & \vdots & \vdots & \vdots & \ddots\\
1 & \xi^{(d_p-1)}_0 & \xi^{(d_p-1)}_1 & (\xi^{(d_p-1)}_0)^2 & \xi^{(d_p-1)}_0 \xi^{(d_p-1)}_1 & (\xi^{(d_p-1)}_1)^2 & \dots
\end{pmatrix}
$$
Note that for any given set of $n$ points $\boldsymbol{\zeta}:=\{\zeta^{(i)}\}_{i=0}^{n-1}$ (which do not have to coincide with the nodal points $\{\xi^{(j)}\}_{j=0}^{d_p-1}$) we can construct the $n\times d_p$ matrix $V(\boldsymbol{\zeta})$ with $V_{ij}(\boldsymbol{\zeta}) = b_j(\zeta^{(i)})$ in the same way. We further define the rank 3 tensor $V^{\partial}(\boldsymbol{\zeta})$ with

$$
V^{\partial}_{ijk}(\boldsymbol{\zeta}):=\frac{\partial b_j}{\partial x_k}(\zeta^{(i)}).
$$
### Tabulation of basis functions
This allows use to *tabulate* the basis functions: for a given set of points $\boldsymbol{\zeta}:=\{\zeta^{(i)}\}_{i=0}^{n-1}$, we have that

$$
T_{ij}(\boldsymbol{\zeta}) := \phi_j(\zeta^{(i)}) = \sum_{m=0}^{d_p-1} c_m^{(j)} b_m(\zeta^{(i)}) = V_{im}(\boldsymbol{\zeta})C_{mj}
$$

or, more compactly:

$$
T(\boldsymbol{\zeta}) = V(\boldsymbol{\zeta}) C
$$

Furthermore, we have for the derivatives

$$
\begin{aligned}
T^\partial_{ijk}(\boldsymbol{\zeta}) &:= \frac{\partial \phi_j}{\partial x_k}(\zeta^{(i)}) 
 = \sum_{m=0}^{d_p-1} c_m^{(j)} \frac{\partial b_m}{\partial x_k}(\zeta^{(i)}) \\
 &= V^\partial_{imk}(\boldsymbol{\zeta})C_{mj}
 \end{aligned}
$$
## Implementation
### Abstract bases class
Since all finite elements share the common functionality that is encapsulated in Ciarlet's definition, we start by writing down an abstract base class, which establishes an interface that all concrete implementations of a finite element need to satisfy. The advantage of this approach is that we do not have to duplicate code that can be shared between all finite element implementation. More specifically, each finite element should provide the following functionality:

* Return the number of nodes associated with each topological entity. For this, we define abstract properties `ndof_per_vertex`, `ndof_per_facet` and `ndof_per_interior` for $\nu_{\text{vertex}}$, $\nu_{\text{facet}}$ and $\nu_{\text{interior}}$ respectively. The base class also contains a property `ndof` which returns $3(\nu_{\text{vertex}}+\nu_{\text{facet}})+\nu_{\text{interior}}$.
* Tabulate the evaluation of all dofs for a given function $\hat{f}$, i.e. compute the vector $(\ell_0(\hat{f}),\ell_1(\hat{f}),\dots,\ell_{d-1}(\hat{f}))^\top$. This is done with the abstract method `tabulate_dofs(fhat)` which gets passed a Python function `fhat`.
* Tabulate the basis functions for a given set of points $\boldsymbol{\zeta}=\{\zeta^{(i)}\}_{i=0}^{n-1}$. This computes the $n\times d$ matrix $T$ with $T_{ij}=\phi_j(\zeta^{(i)})$. This is done with the abstract method `tabulate(zeta)`. If only a single point $\zeta$ is passed to the subroutine it should return a vector of length $d$.
* Tabulate the gradients of all basis functions for a given set of points $\boldsymbol{\zeta}=\{\zeta^{(i)}\}_{i=0}^{n-1}$. This computes the rank 3 tensor $T^\partial$ of shape $n\times d\times 2$ with $T^\partial_{ijk}=\frac{\partial\phi_j}{\partial x_k}(\zeta^{(i)})$. This is done with the abstract method `tabulate_gradient(zeta)`. If only a single point $\zeta$ is passed to the subroutine it should return a matrix of shape $d\times 2$.
* Implement the element dof-map $\mu_{\text{dof}}(E,i,k)$ and its inverse. This is done with the methods `dofmap(entity_type,i,k)` and its inverse `inverse_dofmap(j)`. Since these methods will be called frequently with the same arguments, a [`@functools.cache`](https://docs.python.org/3/library/functools.html#functools.cache) decorator is added to automatically remember  previously used values.

### Concrete implementations
Any concrete implementations of finite elements are obtained by subclassing the `FiniteElement` base class. These concrete classes have to provide concrete implementations of the following methods/properties:
* `ndof_per_vertex`, `ndof_per_facet` and `ndof_per_interior`
* `tabulate_dofs(fhat)` to evaluate the degrees of freedom for a given function
* `tabulate(zeta)` to tabulate the values of the basis functions at a given set of points
* `tabulate_gradient(zeta)` to tabulate the gradients of the basis functions for a given set of points

### Linear element
The bi-linear element is implemented in `LinearElement`

## Exercises
Implement the cubic Lagrange element `CubicElement` ($p=3$) by subclassing the abstract base class `FiniteElement`. The Lagrange points are in this order (see also figure above):

$$
\{\xi^{(j)}\}_{j=0}^{9}=
\left\{
\underbrace{
\begin{pmatrix}0\\[1ex]0\end{pmatrix},
\begin{pmatrix}0\\[1ex]1\end{pmatrix},
\begin{pmatrix}1\\[1ex]1\end{pmatrix}}_{\text{vertices}},
\underbrace{\begin{pmatrix}\frac{2}{3}\\[1ex]\frac{1}{3}\end{pmatrix},
\begin{pmatrix}\frac{1}{3}\\[1ex]\frac{2}{3}\end{pmatrix},
\begin{pmatrix}0\\[1ex]\frac{2}{3}\end{pmatrix},
\begin{pmatrix}0\\[1ex]\frac{1}{3}\end{pmatrix},
\begin{pmatrix}\frac{1}{3}\\[1ex]0\end{pmatrix},
\begin{pmatrix}\frac{2}{3}\\[1ex]0\end{pmatrix}}_{\text{facets}},
\underbrace{\begin{pmatrix}\frac{1}{3}\\[1ex]\frac{1}{3}\end{pmatrix}}_{\text{interior}}\right\}
$$
* You class should store the Lagrange points in an attribute `_nodal_points`
* Your class should contain a method `vandermonde_matrix(zeta,grad=False)` which accepts as an argument a $n\times 2$ matrix of $n$ two-dimensional points. The method should compute the $n\times n_{\text{dof}}$ matrix $V(\boldsymbol{\zeta})$ if `grad=False` and the $n\times n_{\text{dof}}\times 2$ tensor $V^\partial(\boldsymbol{\zeta})$ if `grad=True`.
* Use the `vandermonde_matrix()` method together with `_nodal_points` to construct the coefficient matrix `C`
* Use the coefficient matrix `C` and the `vandermonde_matrix()` method to tabulate the basis functions and their gradients by using the expressions above. You might find the [`numpy.einsum`](https://numpy.org/doc/2.2/reference/generated/numpy.einsum.html) method useful to compute $T^\partial(\boldsymbol{\zeta})$
* Develop a suite of suitable tests to check that your implementation is correct.

### General Lagrange element
The general polynomial element is implemented in the class `PolynomialElement`