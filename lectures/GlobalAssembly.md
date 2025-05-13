# Global assembly

## Interpolation
## Assembly of RHS vector
We compute

$$
\begin{aligned}
b^{(h)}_{i_{\text{global}}} &= \int_\Omega f(x) \Phi^{(h)}_{i_{\text{global}}}(x)\;dx\\
&= \sum_{K\in \Omega_h} \int_K f(x) \Phi_{i_{\text{global}}}(x) \; dx\\
&= \sum_{K\in \Omega_h} \int_K \sum_{i}\widehat{f}_K(\widehat{x}) \phi_i(\widehat{x})\;|\det{J}(\widehat{x})|\;d\widehat{x}\\
&\approx \sum_{K\in \Omega_h} \sum_{i}\sum_q w_q \widehat{f}_K(\xi^{(q)}) \phi_i(\xi^{(q)})\;|\det{J}(\xi^{(q)})|\\
&= \sum_{K\in \Omega_h} \sum_{i}\sum_q w_q \widehat{f}_K(\xi^{(q)}) T_{qi}\;|\det{J}(\xi^{(q)})|
\end{aligned}
$$
with $\widehat{f}_K(\widehat{x}) := f(x)$ and thus
$$
\widehat{f}_K(\xi^{(q)}) = f(x_K^{(q)})
$$
with
$$
x_K^{(q)} = X_K(\xi^{(q)})= \sum_j T^\times_{qj} X_{j_{\text{global}}(K,j)}
$$
The Jacobian is given by
$$
J_{ab}(\xi^{(q)}) = \frac{\partial (X_K)_a }{\partial x_b}(\xi^{(q)})
= \sum_j X_{j_{\text{global}}} \frac{\partial (\phi^\times_j)_a }{\partial x_b}(\xi^{(q)})
= \sum_j X_{j_{\text{global}}(K,j)} T^{\times\partial}_{qjab}
$$

### Algorithm
1. Initialise $\boldsymbol{b}^{(h)} \gets \boldsymbol{0}$
1. For all cells $K$ **do**:
1. $~~~~$ Extract the coordinate dof-vector $\overline{X}$ with $\overline{X}_j = X_{j_\text{global}(K,j)}$
1. $~~~~$ For all quadrature points $q$ **do**:
1. $~~~~~~~~$ Compute the determinant $D_q$ of the Jacobian $J(\xi^{(q)})$ with $J_{ab}(\xi^{(q)}) = \sum_j \overline{X}_j T^{\times\partial}_{qjab}$
1. $~~~~~~~~$ Compute $x_K^{(q)} = \sum_j T^\times_{qj} \overline{X}_j$ and evaluate $F_q = f(x_K^{(q)})$
1. $~~~~$ **end do**
1. $~~~~$ For all local dof-indices $i$ **do**:
1. $~~~~~~~~$ Increment $b_{i_{\text{global}}}\gets b_{i_{\text{global}}} + \sum_q w_q F_q T_{qi} D_q$
1. $~~~~$ **end do**
1. **end do**


## Assembly of LHS matrix
We compute

$$
\begin{aligned}
A^{(h)}_{i_{\text{global}},j_{\text{global}}} &= \int_\Omega \left(\kappa \nabla \Phi^{(h)}_{i_{\text{global}}}(x) \cdot\nabla\Phi^{(h)}_{j_{\text{global}}}(x)+\omega \Phi^{(h)}_{i_{\text{global}}}(x)\Phi^{(h)}_{j_{\text{global}}}(x)\right)dx\\
&= \sum_{K\in\Omega_h}\int_K \left(\kappa \nabla \Phi^{(h)}_{i_{\text{global}}}(x) \cdot\nabla\Phi^{(h)}_{j_{\text{global}}}(x)+\omega \Phi^{(h)}_{i_{\text{global}}}(x)\Phi^{(h)}_{j_{\text{global}}}(x)\right)dx\\
&= \sum_{K\in \Omega_h}\int_K \sum_{ij} \left(\kappa J^{-\top}(\widehat{x}) \widehat{\nabla} \phi_i(\widehat{x})\cdot J^{-\top}(\widehat{x})\phi_j(\widehat{x}) + \omega\phi_i(\widehat{x})\phi_j(\widehat{x})\right)|\det{J}(\widehat{x})|d\widehat{x}\\
&\approx \sum_{K\in \Omega_h}\int_K \sum_{ij} w_q \left(\kappa \widehat{\nabla} \phi_i(\xi^{(q)})(J^{\top}(\xi^{(q)}) J(\xi^{(q)}))^{-1}\phi_j(\xi^{(q)}) + \omega\phi_i(\xi^{(q)})\phi_j(\xi^{(q)})\right)|\det{J}(\xi^{(q)})|d\widehat{x} \\
&= \sum_{K\in \Omega_h} \sum_{ij}\sum_q w_q  \left(\kappa \sum_{k\ell}T_{qik}(J^{(2)}_q)_{k\ell} T_{qj\ell} +\omega T_{qi}T_{qj}\right)|\det{J}(\xi^{(q)})|
\end{aligned}
$$

$$
J^{(2)}_{q} =  \left(J^{\top}(\xi^{(q)}) J(\xi^{(q)})\right)^{-1}
$$

1. Initialise $A^{(h)} \gets 0$
1. For all cells $K$ **do**:
1. $~~~~$ For all quadrature points $q$ **do**:
1. $~~~~$ Extract the coordinate dof-vector $\overline{X}$ with $\overline{X}_j = X_{j_\text{global}(K,j)}$
1. $~~~~~~~~$ Compute the Jacobian $J(\xi^{(q)})$ with $J_{ab}(\xi^{(q)}) = \sum_j \overline{X}_j T^{\times\partial}_{qjab}$
1. $~~~~~~~~$ Compute the determinant $D_q$ of $J(\xi^{(q)})$
1. $~~~~~~~~$ Compute the matrix $J^{(2)}_{q} =  \left(J^{\top}(\xi^{(q)}) J(\xi^{(q)})\right)^{-1}$
2. $~~~~$ **end do**
3. $~~~~$ For all local dof-indices $i$ **do**:
4. $~~~~~~~~$ For all local dof-indices $j$ **do**:
5. $~~~~~~~~~~~~$ Increment $A^{(h)}_{i_{\text{global}},j_{\text{global}}}\gets A^{(h)}_{i_{\text{global}},j_{\text{global}}} + w_q  \left(\kappa \sum_{k\ell}T_{qik}(J^{(2)}_q)_{k\ell} T_{qj\ell} +\omega T_{qi}T_{qj}\right)D_q$
6. $~~~~~~~~$ **end do**
7. $~~~~$ **end do**
8.  **end do**
