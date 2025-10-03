----

*&#169; Eike Mueller, University of Bath 2025. These notes are copyright of Eike Mueller, University of Bath. They are provided exclusively for educational purposes at the University and are to be downloaded or copied for your private study only. Further distribution, e.g. by upload to external repositories, is prohibited. html generated with [pandoc](https://pandoc.org/) using [easy-pandoc-templates](https://github.com/ryangrose/easy-pandoc-templates) under the [GPL-3.0.1 license](https://github.com/ryangrose/easy-pandoc-templates?tab=GPL-3.0-1-ov-file#readme)*

----
## Implementation

The code can be made very compact by implementing tensor contractions with the [`np.einsum()` method](https://numpy.org/doc/stable/reference/generated/numpy.einsum.html). Consider, for example, the assembly of the left-hand side siffness matrix $A^{(h)}_{\ell k}$ in `assemble_lhs()`. From the lectures we know that

$$
\begin{aligned}
A^{(h)}_{\ell k} &=  \kappa \sum_{q=0}^{N_q-1}\sum_{a=0}^{d-1} w_q  T^\partial_{q\ell a} T^\partial_{qka} 
+\omega \sum_{q=0}^{N_q-1} w_qT_{q\ell} T_{qk}
\end{aligned}
$$

with the quadrature weights $w_q$ and the tensors $T$, $T^\partial$ of shapes $(N_q,\nu)$, $(N_q,\nu,d)$ (with $d=2$) respectively where $T_{q\ell} = \phi_\ell(\zeta^{(q)})$ and $T_{q\ell a}^{\partial} = \frac{\partial\phi_\ell}{\partial x_a}(\zeta^{(q)})$.

In the code, we start by constructing a quadrature rule and construct the vector $\boldsymbol{w}\in\mathbb{R}^{N_q}$ which holds the weights:
```Python
quad = GaussLegendreQuadratureReferenceTriangle(n_q)
w_q = quad.weights
```
The quadrature points $\zeta$ are stored in an array `w_q` of shape $(N_q,2)$, which can be passed to the `tabulate()` and `tabulate_gradient()` methods of the finite element

```Python
zeta_q = np.asarray(quad.nodes)
phi = element.tabulate(zeta_q)
grad_phi = element.tabulate_gradient(zeta_q)
```

As a result, we obtain an array `phi` of shape $(N_q,\nu)$ which stores the rank-2 tensor $T$ and an array `grad_phi` of shape $(N_q,\nu,d)$ which stores the rank-3 tensor $T^\partial$. The stiffness matrix $A^{(h)}_{\ell k}$ is then constructed with two suitable tensor contractions involving `phi`, `grad_phi` and `w_q`:

```Python
stiffness_matrix = kappa * np.einsum("q,qik,qjk->ij", w_q, grad_phi, grad_phi) 
                 + omega * np.einsum("q,qi,qj->ij", w_q, phi, phi)
```

The complete code for assembling the stiffness matrix can be found in [algorithms.py](algorithms.py). The main program has been implemented in [driver.py](driver.py).

## Numerical experiments

### Error norm
The following table shows the $L_2$ norm of the error for polynomial degrees $p=1,2,3,4$

| degree $p$ | error norm $\|\|e_h\|\|_{L_2(\widehat{K})}$ |
| ---------- | ---------------------------- |
|         1  |          $9.06\cdot 10^{-2}$ |
|         2  |          $2.75\cdot 10^{-2}$ |
|         3  |          $8.38\cdot 10^{-3}$ |
|         4  |          $2.82\cdot 10^{-3}$ |

### Plot of solution and error
The following plots visualise the solution and error for polynomial degrees $p=1,2,3,4$

#### Solution and error for $p=1$
![Solution and error for $p=1$](triangle_solution_01.png)

#### Solution and error for $p=2$
![Solution and error for $p=2$](triangle_solution_02.png)

#### Solution and error for $p=3$
![Solution and error for $p=3$](triangle_solution_03.png)

#### Solution and error for $p=4$
![Solution and error for $p=4$](triangle_solution_04.png)