# Unstructured meshes
In general, we might want to solve PDEs on a $d$-dimensional manifold $\Omega\subset \mathbb{R}^D$ that is embedded in $D$ dimensional space. For example, we might want to solve the Navier-Stokes equations on the surface of a sphere. The manifold is then approximated by a mesh, which can be described as a collection of topological entities. For example, if $d=2$, the mesh will consist of zero-dimensional vertices, one-dimensional edges and two-dimensional triangular cells (although it is also possible to use more general polygonal cells we do not consider this here). In general, the co-dimension $c$ of a $d'$-dimensional mesh entity is given by $c=d-d'$. In the following we will only consider the case $d=D=2$. In this case we have:

| topological entity  | dimension $d'$ | co-dimension $c$ |
| ------------------- | -------------- | ---------------- |
| cell (triangle) $K$ | $2$            | $0$              |
| facet (edge) $F$    | $1$            | $1$              |
| vertex $v$          | $0$            | $2$              |

The following figure shows a two-dimensional mesh in which all topological entities are labelled by their co-dimensional and a unique number.

![simple triangular mesh](figures/simple_mesh.svg)

## Topology
The mesh topology is defined by which entities are connected to which other entities. For example, each cell has exactly three facets and each facet is defined by exactly two vertices.

This information can be encoded in two matricies: An $n_{\text{cell}}\times 3$ matrix $I^{F\gets K}$ with

$$
I^{F\gets K}_{ij} = \text{index of $j$-th facet of cell $i$}
$$

and an $n_{\text{facet}}\times 2$ matrix $I^{v\gets F}$ with

$$
I^{v\gets F}_{jk} = \text{index of $k$-th vertex of facet $j$}.
$$

For convenience, we can also use $I^{F\gets K}$ and $I^{v\gets F}$ to construct the $n_{\text{cell}}\times 3$ matrix $I^{v\gets K}$ with

$$
I^{v\gets K}_{ik} = \text{index of $k$-th vertex of cell $i$}
$$

More explicitly, the $i$-th column of this matrix is given by

$$
\begin{aligned}
I^{v\gets K}_{i,0} &= I_{j_0,0}^{v\gets F}\quad\text{with}\;j_0=I_{i,2}^{F\gets K}\\
I^{v\gets K}_{i,1} &= I_{j_1,0}^{v\gets F}\quad\text{with}\;j_1=I_{i,0}^{F\gets K}\\
I^{v\gets K}_{i,2} &= I_{j_2,0}^{v\gets F}\quad\text{with}\;j_2=I_{i,1}^{F\gets K}.
\end{aligned}
$$

For example, the matrices $I^{F\gets K}$, $I^{v\gets F}$ and $I^{v\gets K}$ for the simple mesh shown above are given by

$$
\begin{aligned}
I^{F\gets K} &= \begin{pmatrix}
0 & 5 & 1 & 0 & 3 \\
1 & 3 & 7 & 5 & 4 \\
8 & 2 & 2 & 9 & 6
\end{pmatrix}^{\top}\\
I^{v\gets F} &= \begin{pmatrix}

0 & 1 & 2 & 5 & 2 & 1 & 3 & 4 & 4 & 5 \\
1 & 4 & 1 & 2 & 3 & 5 & 5 & 2 & 0 & 0
\end{pmatrix}^{\top}\\
I^{F\gets K} &= \begin{pmatrix}
4 & 2 & 2 & 5 & 3 \\
0 & 1 & 1 & 0 & 5 \\
1 & 5 & 4 & 1 & 2
 \\

\end{pmatrix}^{\top}

\end{aligned}
$$

## Grid cells and reference elements
## Encoding geometry information
## Implementation in Python