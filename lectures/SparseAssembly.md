# Assembly of stiffness matrix

### Algorithm: create sparsity structure
1. Set $\mathcal{J}_{\ell_{\text{global}}} = \emptyset$ for all $\ell_{\text{global}}=0,1,\dots,n_{\text{dof,global}}-1$
2. **for all** cells $K$ with index $i$
3. $~~~~$ Set $\mathcal{L}^{(i)} = \{\ell_{\text{global}}(i,\ell)\;\text{for}\;\ell=0,1,\dots,n_{\text{dof,local}}-1\}$, the set of global indices of dofs associated with cell $K$
4. $~~~~$ **for all** $\ell_{\text{global}}\in \mathcal{L}^{(i)}_{\text{global}}$ **do** 
5. $~~~~~~~~$ Update $\mathcal{J}_{\ell_{\text{global}}} \gets \mathcal{J}_{\ell_{\text{global}}} \cup \mathcal{L}_{\text{global}}^{(i)}$
6. $~~~~$ **end do** 
7. **end do** 
8. Initialise $R = [0,0,\dots,0]\in \mathbb{R}^{n+1}$
9. Initialise $J = []$
10. **for** $\ell_{\text{global}}=0,1,\dots,n_{\text{dof,global}}-1$ **do**
11. $~~~~$ Append $\mathcal{J}_{\ell_{\text{global}}}$ to $J$
12. $~~~~$ Set $R_{\ell_{\text{global}}+1} = R_{\ell_{\text{global}}} +\left|\mathcal{J}_{\ell_{\text{global}}}\right|$ 
13. **end do**

    