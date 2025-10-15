----

*&#169; Eike Mueller, University of Bath 2025. These notes are copyright of Eike Mueller, University of Bath. They are provided exclusively for educational purposes at the University and are to be downloaded or copied for your private study only. Further distribution, e.g. by upload to external repositories, is prohibited. html generated with [pandoc](https://pandoc.org/) using [easy-pandoc-templates](https://github.com/ryangrose/easy-pandoc-templates) under the [GPL-3.0.1 license](https://github.com/ryangrose/easy-pandoc-templates?tab=GPL-3.0-1-ov-file#readme)*

----
# Backsubstitution

## Pseudocode
The solution vector $\boldsymbol{u}$ of the system $A\boldsymbol{u}=\boldsymbol{b}$ can be computed as follows if $A$ is an upper triangular matrix:

1. **for** $i=n-1,n-2,\dots,0$ **do**
2. $~~~~$ Set $r = b_i$
3. $~~~~$ **for** $j=i+1,i+2,\dots,n-1$ **do**
4. $~~~~~~~~$ Update $r\gets r - A_{ij}\cdot u_j$
5. $~~~~$ **end for**
6. $~~~~$ Set $u_i = r/A_{ii}$
7. **end for**

## Computational complexity
At the $i$-step of the algorithm we compute

$$
u_i = \left(b_i-\sum_{j=i+1}^{n-1}\cdot A_{ij} u_j \right)/A_{ii}
$$

This requires $n-1-i$ multiplications of the form $A_{ij}\cdot u_j$ as well as $n-1-i$ subtractions/additions (one for each term in the sum) and one division by $A_{ii}$. We need to do this for $i=n-1,n-1,\dots,0$ (i.e. each iteration of the outer loop) and hence the total number of operations is

$$
\begin{aligned}
C_{\text{backsub}}(n) &= \sum_{i=0}^{n-1} (2(n-1-i)+1) \\
&= \sum_{j=0}^{n-1} (2j+1)\\
&= n(n-1)+n = n^2 = \mathcal{O}(n^2)
\end{aligned}
$$

# Solving triangular systems

## Reduction to upper tridiagonal system
Consider the $n\times n$ tridiagonal matrix $A$. To remove all entries below the diagonal in the first column, we need to scale the first row by $\rho = -A_{10}/A_{00}$ and add it to the second row. We also need to replace $b_1 \mapsto b_1 - A_{10}/A_{00}\cdot b_0 = b_1+\rho\cdot b_0$. Since $A_{01}$ is the only non-zero entry in the first row which is not on the diagonal, this requires the following operations:

* compute $\rho = -A_{10}/A_{00}$ $\Rightarrow$ **$\boldsymbol{1}$ division**
* Replace $A_{11} \mapsto A_{11} + \rho\cdot A_{01}$ $\Rightarrow$ **$\boldsymbol{1}$ multiplication** and **$\boldsymbol{1}$ addition**, since we can ignore the first entry in the second row: this entru will be set to zero by construction
* update $b_1 \gets b_1 + \rho\cdot b_1$ $\Rightarrow$ **$\boldsymbol{1}$ multiplication** and **$\boldsymbol{1}$ addition**

Next, we apply this process to the tridiagonal $(n-1)\times (n-1)$ submatrix in the lower right corner. Continuing recursively $n-1$ times until the remaining matrix is of size $1\times 1$, we obtain an upper triangular matrix with exactly one superdiagonal. Since each of the $n-1$ steps requires 5 floating point operations, the total number of operations is therefore

$$
n_{\text{reduction}}(n) = 5(n-1) = 5n-5
$$

## Backsubstitution
To solve the resulting system with backsubstitution, we can use the above algorithm, but only need to consider the case $j=i+1$ in the innermost loop. This results in:

1. **for** $i=n-1,n-2,\dots,0$ **do**
2. $~~~~$ **if** $i=n-1$ **then**
3. $~~~~~~~~$ Set $u_i = b_i/A_{ii}$
4. $~~~~$ **else**
5. $~~~~~~~~$ Set $r = b_i - A_{i,i+1}\cdot u_{i+1}$
6. $~~~~~~~~$ Set $u_i = r/A_{ii}$
7. $~~~~$ **end if**
7. **end for**

This requires
$$
n_{\text{backsub}}(n) = 3(n-1) + 1 = 3n - 2
$$
operations. The total number of operations required to solve a triangular system is therefore

$$
n_{\text{solve}}(n) = 8(n-1) + 1 = 8n - 7 = \mathcal{O}(n).
$$

This should be compared to the $\mathcal{O}(n^3)$ computational complexity which is required to solve a general linear system of the form $A\boldsymbol{u}=\boldsymbol{b}$.