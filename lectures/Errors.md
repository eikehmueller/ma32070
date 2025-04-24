# Sources of errors

## Floating point numbers

## Rounding errors

### Simple example

Consider the linear system

$$
\begin{pmatrix}
0 & 1 \\ 1 & 1 
\end{pmatrix}
\begin{pmatrix}
x_0 \\ x_1
\end{pmatrix}
=
\begin{pmatrix}
1 \\ 0
\end{pmatrix}
$$
The solution is $x_0=-1$, $x_1=+1$. Now consider a small perturbation of this problem, namely
$$
\begin{pmatrix}
-10^{-20} & 1 \\ 1 & 1 
\end{pmatrix}
\begin{pmatrix}
x_0 \\ x_1
\end{pmatrix}
=
\begin{pmatrix}
1 \\ 0
\end{pmatrix}
$$
It is easy to see that the solution of the perturbed problem is $x_0=-\frac{1}{1+10^{-20}}$, $x_1=+\frac{1}{1+10^{-20}}$, which is very close to the solution of the unperturbed system.

Let us solve the perturbed system numerically. For this, we write it as

$$
\begin{aligned}
-10^{-20} x_0 + x_1 &= 1\\
x_0 + x_1 &= 0
\end{aligned}
$$
To eliminate $x_0$ from the second equation, we multiply the first equation with $10^{20}$ and add it to the second equation to obtain
$$
\begin{aligned}
-10^{-20} x_0 + x_1 &= 1\\
(1+10^{20}) x_1 &= 10^{20}
\end{aligned}
$$
Now, since $10^{20}\gg 1$, we can replace $1+10^{20}$ by $10^{20}$ to obtain
$$
\begin{aligned}
-10^{-20} x_0 + x_1 &= 1\\
10^{20} x_1 &= 10^{20},
\end{aligned}
$$
which immediately implies $x_1 = 1$. Inserting this into the first equation gives
$$
\begin{aligned}
-10^{-20} x_0 + 1 &= 1
\end{aligned}
$$
and therefore $x_0=0$. Altogether we find $x_0=0$, $x_1=1$. This is very different from the exact solution $x_0=-\frac{1}{1+10^{-20}}$, $x_1=+\frac{1}{1+10^{-20}}$! Hence, although the rounding we performed in the numerical solution procedure appears to be innocent, we get a completely wrong solution. In this case, this problem can be fixed by using a slightly different solution procedure: subtract the first equation from the second equation to obtain $(1+10^{-20}) x_0 = -1$, which can be safely approximated by $x_0=-1$. Then use this in the second equation $x_0+x_1=0$ to conclude $x_1=+1$. Now this is very close to the exact solution of the perturbed system. It turns out that solving the perturbed linear system with [`np.linalg.solve()`](https://numpy.org/doc/2.2/reference/generated/numpy.linalg.solve.html) will also give a good approximation. Unfortunately, this is not always the case, as the following simple example demonstrates.

### More complicated example

For a given positive $\epsilon>0$ consider the following $2\times 2$ symmetric matrix

$$
\begin{aligned}
A &= \begin{pmatrix}
1+\frac{1}{\sqrt{2}} + \frac{2+\sqrt{2}}{4}\epsilon & \frac{1}{\sqrt{2}}-\frac{\sqrt{2}}{4}\epsilon\\[1ex]
\frac{1}{\sqrt{2}}-\frac{\sqrt{2}}{4}\epsilon & 1-\frac{1}{\sqrt{2}}+\frac{2-\sqrt{2}}{4}\epsilon
\end{pmatrix}
\\
&\approx
\begin{pmatrix}
1.7071067811865475 + 0.8535533905932737 \cdot\epsilon & 0.7071067811865475 - 0.3535533905932738 \cdot\epsilon \\
0.7071067811865475 - 0.3535533905932738 \cdot\epsilon & 0.2928932188134525 + 0.1464466094067262\cdot\epsilon
\end{pmatrix}
\end{aligned}
$$

and vector

$$
\boldsymbol{b} = \begin{pmatrix}
\frac{\sqrt{2+\sqrt{2}}+\sqrt{2-\sqrt{2}}}{2}\\[1ex]
\frac{\sqrt{2+\sqrt{2}}-\sqrt{2-\sqrt{2}}}{2}
\end{pmatrix}
\approx
\begin{pmatrix}
1.3065629648763766\\0.5411961001461970
\end{pmatrix}
$$

It can be shown that independent of $\epsilon$ the exact solution of the linear system $A\boldsymbol{u}=\boldsymbol{b}$ is given by

$$
\boldsymbol{u}_{\text{exact}} = 
\begin{pmatrix}
\frac{\sqrt{2 - \sqrt{2}}}{2}\\[1ex]
\frac{\sqrt{2 + \sqrt{2}}}{2}
\end{pmatrix}
\approx 
\begin{pmatrix}
0.3826834323650897\\0.9238795325112867
\end{pmatrix}
$$

Hence, we would expect that if we solve the linear system with a numerical method, we would get a solution that is at least close to the exact solution. The following table shows the solution $\boldsymbol{u}$ of $A\boldsymbol{u}=\boldsymbol{b}$ that is obtained with
```
u = np.linalg.solve(A,b)
```
for different values of $\epsilon$. The final column shows the relative error $\|\boldsymbol{u}-\boldsymbol{u}_{\text{exact}}\|_2|/\|\boldsymbol{u}_{\text{exact}}\|_2$

| $\epsilon$        | solution $\boldsymbol{u}$                                                | relative error $\|\|\boldsymbol{u}-\boldsymbol{u}_{\text{exact}}\|\|_2/\|\|\boldsymbol{u}_{\text{exact}}\|\|_2$ | condition number $\kappa$ |
| ----------------- | ------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------- | ------------------------- |
| $10^{-3}$         | $\begin{pmatrix}0.3826834323650560 \\ 0.9238795325113685\end{pmatrix}$   | $8.8408\cdot 10^{-14}$                                                                                          | $4\cdot 10^3$             |
| $10^{-6}$         | $\begin{pmatrix}0.3826834323171860 \\0.9238795326269368\end{pmatrix}$    | $1.2518\cdot 10^{-10}$                                                                                          | $4\cdot 10^6$             |
| $10^{-9}$         | $\begin{pmatrix}0.3826834290672392 \\0.9238795404730024\end{pmatrix}$    | $8.6177\cdot 10^{-9}$                                                                                           | $4\cdot 10^9$             |
| $10^{-12}$        | $\begin{pmatrix}0.3826776593455781 \\ 0.9238934698132879\end{pmatrix}$   | $1.5086\cdot10^{-5}$                                                                                            | $4\cdot 10^{12}$          |  |
| $10^{-15}$        | $\begin{pmatrix}0.3888090807546386 \\0.9090909090909092\end{pmatrix}$    | $1.6007\cdot10^{-2}$                                                                                            | $4\cdot 10^{15}$          |
| $8\cdot 10^{-16}$ | $\begin{pmatrix} 0.2919799363037854 \\1.1428571428571428\end{pmatrix}$   | $2.3702\cdot 10^{-1}$                                                                                           | $5\cdot 10^{15}$          |
| $4\cdot 10^{-16}$ | $\begin{pmatrix} -0.0630602600160104 \\ 2.0000000000000000\end{pmatrix}$ | $1.1648$                                                                                                        | $10^{16}$                 |

Clearly, the errors increase for smaller values of $\epsilon$. This is related to the fact that the condition number of the matrix increases: The eigenvalues of $A$ are

$$
\lambda_{\pm} = 1+\frac{\epsilon}{2}\pm\sqrt{1-\frac{\epsilon^2}{4}} \approx \{2,\frac{\epsilon}{2}\} \qquad\text{for $\epsilon\ll 1$}
$$

and hence the condition number, which is the ratio of the largest and smallest eigenvalue, is

$$
\kappa = \frac{\lambda_+}{\lambda_-} = \frac{1+\frac{\epsilon}{2}+\sqrt{1-\frac{\epsilon^2}{4}}}{1+\frac{\epsilon}{2}-\sqrt{1-\frac{\epsilon^2}{4}}} \approx \frac{4}{\epsilon}
$$

## Backward error analysis