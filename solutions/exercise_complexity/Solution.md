<div align="center">
  <p style="font-size:32px;">MA32070 Exercise 4: Model solution</p>
</div>

----

*&#169; Eike Mueller, University of Bath 2025. These notes are copyright of Eike Mueller, University of Bath. They are provided exclusively for educational purposes at the University and are to be downloaded or copied for your private study only. Further distribution, e.g. by upload to external repositories, is prohibited. html generated with [pandoc](https://pandoc.org/) using [easy-pandoc-templates](https://github.com/ryangrose/easy-pandoc-templates) under the [GPL-3.0.1 license](https://github.com/ryangrose/easy-pandoc-templates?tab=GPL-3.0-1-ov-file#readme)*

----
## Pseudocode for backsubstitution
The solution vector $\boldsymbol{u}$ of the system $A\boldsymbol{u}=\boldsymbol{b}$ can be computed as follows if $A$ is an upper triangular matrix:

1. **for** $i=n-1,n-2,\dots,0$ **do**
2. $~~~~$ Set $r = b_i$
3. $~~~~$ **for** $j=i+1,i+2,\dots,n-1$ **do**
4. $~~~~~~~~$ Update $r\gets r - A_{ij}u_j$
5. $~~~~$ **end for**
6. $~~~~$ Set $u_i = r/A_{ii}$
7. **end for**

## Computational complexity
At the $i$-step of the algorithm we compute

$$
u_i = \left(b_i-\sum_{j=i+1}^{n-1} A_{ij} u_j \right)/A_{ii}
$$

This requires $n-1-i$ multiplications $A_{ij} u_j$, $n-1-i$ subtractions/additions and one division. We need to do this for $i=n-1,n-2,\dots,0$, hence the total number of operations is

$$
\begin{aligned}
C_{\text{backsub}}(n) &= \sum_{i=0}^{n-1} (2(n-1-i)+1) \\
&= \sum_{j=0}^{n-1} (2j+1)\\
&= n(n-1)+n = n^2
\end{aligned}
$$

## Solving triangular systems