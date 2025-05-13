# Function spaces

## Grid cells and reference elements
We can now construct a function space $V_h\subset H^1(\Omega_h)$ on the domain $\Omega_h$ defined by the mesh as follows. First, we assume that each cell $K$ of the mesh is the image of the reference cell $\widehat{K}$ under a map $X_K$. Next, assume that the function $u|_K$ in each cell is defined such that its pullback under the map $X_K$ is a polynomial. More specifically we set
$$
V_h := \{u\in H^1(\Omega_h): u|_K(x) = \widehat{u}_K(\widehat{x})\quad \text{for some $\widehat{u}_K\in\mathcal{P}_p(\widehat{K})$ with $x=X_K(\widehat{x})$ for all cells $K\in\Omega_h$}\}
$$

![pullback to reference element](figures/reference_mapping.svg)

We can now use the basis for $\mathcal{P}_p(\widehat{K})$ that we constructed in one of the previous lectures to represent the functions $\widehat{u}_K$. However, care has to be taken to ensure that the function is continuous across facets and at the vertices. To guarantee this, we need to think carefully about the arrangement of unknowns on the mesh.

### Arrangement of unknowns
Assume that we have a finite element with $\nu_{\text{vertex}}$ unknowns per vertex, $\nu_{\text{facet}}$ unknowns per facet and $\nu_{\text{interior}}$ unknowns per cell. 
Let $N_{\text{vertex}} = n_{\text{vertex}}\cdot \nu_{\text{vertex}}$, $N_{\text{facet}} = n_{\text{facet}}\cdot \nu_{\text{facet}}$ and $N_{\text{interior}} = n_{\text{cell}}\cdot \nu_{\text{interior}}$ be the total number of unknowns associated with vertices, facets and cell interiors respectively.

We number the unknowns by using the first $N_{\text{vertex}}$ indices for unknowns associated with vertices, the next $N_{\text{facet}}$ indices for unknowns associated with cells and the remaining $N_{\text{interior}}$ indices for unknowns associated with cell interiors. More specifically:

* unknowns with indices $k\cdot \nu_{\text{vertex}},\dots,(k+1)\cdot \nu_{\text{vertex}}-1$ are associated with the $k$-th vertex
* unknowns with indices $N_{\text{vertex}}+j\cdot \nu_{\text{facet}},\dots,N_{\text{vertex}}+(j+1)\cdot \nu_{\text{facet}}-1$ are associated with the $j$-th facet
* unknowns with indices $N_{\text{vertex}}+N_{\text{facet}}+i\cdot \nu_{\text{interior}},\dots,N_{\text{vertex}}+N_{\text{facet}}+(i+1)\cdot \nu_{\text{interior}}-1$ are associated with the $i$-th cell

The facet-unknowns are ordered along the orientation of each edge.

The following figure shows an example for the $p=3$ Lagrange element. For this mesh there are $N_{\text{vertex}}=6$ unknowns associated with the vertices, $N_{\text{facet}}=10\cdot 2=20$ unknowns associated with the facets and $N_{\text{interior}} = 5$ unknowns associated with the cell interiors.

![simple triangular mesh with unknowns](figures/simple_mesh_with_dof_numbers.svg)

### Local-to-global map
On each cell with index $i$ we need to map the $\ell$-th local dof-index to the global dof-index $\ell_{\text{global}}$. Note that this map is surjective but not injective since the unknowns associated with the facets and vertices are shared.

The map $(i,\ell) \mapsto \ell_{\text{global}}$ is realised in two steps:

1. For the given $\ell$, we first work out the entity-type $E$ (cell, facet or vertex) it is associated with, as well as the index $j$ of this entity on the reference triangle and the index $k$ of the dof on that reference entity. Recall that for the finite element we have the dof-map $\mu_{\text{dof}}$ that $\ell = \mu_{\text{dof}}(E,j,k)$, so we can obtain $E$, $j$ and $k$ from the inverse of this map.
2. Next, we map this to the global index taking into account the arrangement of unknowns described above:
   1. If $E=\text{vertex}$: $\ell_{\text{global}} = v\cdot \nu_{\text{vertex}}+k$ where $v=I^{v\gets K}_{ij}$ is the global index of the vertex with local index $j$.
   2. If $E=\text{facet}$: $\ell_{\text{global}} = N_{\text{vertex}}+f\cdot \nu_{\text{facet}}+\widetilde{k}$ where $f=I^{F\gets K}_{ij}$ is the global index of the facet with local index $j$. Note that the orientation of the facet does not necessarily match the orientation on the reference element $\widehat{K}$. To take this into account, set $\widetilde{k}=k$ if the orientations agree and set $\widetilde{k} = \nu_{\text{facet}}-1-k$ otherwise.
   3. If $E=\text{cell}$: $\ell_{\text{global}} = N_{\text{vertex}}+N_{\text{facet}}+i\cdot \nu_{\text{interior}}+k$

## Encoding geometry information
The geometry of the mesh is encoded in the functions $X_K$. We can combine the $X_K$ in each cell into a function $X$ which is defined on the entire mesh. The crucial observation is now that each component of $X_K$ can be represented by a function in a finite element space $W_h$ which is defined in the same way as $V_h$ (but possibly with a different polynomial degree). Hence, we can define a coordinate function 

$$
X \in W_h^{\times} = W_h \times W_h
$$

The space $W_h^{\times}$ is a vector function space. Let $w=(w_0,w_1):\widehat{K}\rightarrow \mathbb{R}^2$ be a vector-valued function. Then the degrees of freedom $\ell^\times$ of $W^\times_h$ are given by:

$$
\ell^\times_j(w) = \begin{cases} 
\ell_{j/2}(w_0) & \text{for $j$ even}\\
\ell_{(j-1)/2}(w_1) & \text{for $j$ odd}\\
\end{cases}
$$

The vector-valued basis functions $\phi^\times_j$ are
$$
\phi^\times_j(x) = 
\begin{cases}
\begin{pmatrix}\phi_{j/2}(x)\\0\end{pmatrix} & \text{for $j$ even}\\[2ex]
\begin{pmatrix}0 \\\phi_{(j-1)/2}(x)\end{pmatrix} & \text{for $j$ odd}
\end{cases}
$$

### Tabulation of basis functions
As for scalar-valued function spaces we can *tabulate* the basis functions. For a given set of points $\boldsymbol{\zeta}:=\{\zeta^{(i)}\}_{i=0}^{n-1}$, we obtain the $n \times \nu^\times \times 2$ matrix $T^\times$ with

$$
T^\times_{ij\ell}(\boldsymbol{\zeta}) := (\phi^\times_j(\zeta^{(i)}))_\ell = 
\begin{cases}
\phi_{j/2}(\zeta^{(i)}) = T_{i,j/2}& \text{if $j$ even and $\ell=0$} \\
\phi_{(j-1)/2}(\zeta^{(i)}) = T_{i,(j-1)/2} & \text{if $j$ odd and $\ell=1$} \\
0 & \text{otherwise}
\end{cases}
$$

Furthermore, the derivatives are collected in the $n \times \nu^\times \times 2\times 2$ matrix $T^{\times\partial}$ matrix with

$$
T^{\times\partial}_{ij\ell k}(\boldsymbol{\zeta}) := \frac{(\partial \phi^\times_j)_\ell}{\partial x_k}(\zeta^{(i)})  = \begin{cases}
\frac{\phi_{j/2}}{\partial x_k}(\zeta^{(i)}) = T_{i,j/2,k}& \text{if $j$ even and $\ell=0$} \\
\frac{\partial\phi_{(j-1)/2}}{\partial x_k}(\zeta^{(i)}) = T_{i,(j-1)/2,k} & \text{if $j$ odd and $\ell=1$} \\
0 & \text{otherwise}
\end{cases}
$$

We can now write 

$$
X_K(\xi^{(i)}) = \sum_j X_j \phi_j^\times(\xi^{(i)}) = \sum_j T^\times_{ij} X_j
$$

and the Jacobian is

$$
J_{ab}(\xi^{(i)}) = \frac{\partial (X_K)_a }{\partial x_b}(\xi^{(i)})
= \sum_j X_j \frac{\partial (\phi^\times_j)_a }{\partial x_b}(\xi^{(i)})
= \sum_j X_j T^{\times\partial}_{ijab}
$$


## Implementation in Python