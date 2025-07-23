# Error analysis

## Sources of error
When solving a problem in Scientific Computing, there are several sources of error:

#### Modelling error
When modelling a physical phenomenon, we need to pick a set of equations. For example, we might want to use the Navier-Stokes equations to model fluid flow in the atmosphere. Since most equations are an approximation of the real physics, this will inevitably introduce modelling errors.

#### Discretisation error
To solve the chosen system of equations they need to be discretised so that they can be solved on a computer. The finite element discretisation will introduce errors that are typically of the form $Ch^\alpha$ for some positive constants $C$, $\alpha$ where $h$ is the grid spacing. The error can be reduced by refining the compututational grid or by choosing a better discretisation which might lead to smaller $C$ and larger $\alpha$.

#### Computational (or rounding) error
Since a computer can only perform inexact arithmetic for real numbers, the results will only be accurate up to rounding errors.

Obviously, it is crucial to minimise the total error, which is made up of the three components above. Modelling errors are discussed elsewhere and beyond the scope of this course, in which we will concentrate on the PDE $-\kappa \Delta u + \omega u = f$. A detailled analysis of finite element discretisation errors is presented in the course on "Numerical solution of elliptic PDEs". In the following we focus on rounding errors.

## Results from numerical experiment
As a motivation, consider the solution of our model equation $-\kappa \Delta u + \omega u = f$ on the reference triangle for $\kappa = 0.9$, $\omega = 0.4$. The boundary conditions and right-hand side were chosen such that the exact solution is given by $u_{\text{exact}}(x) = \exp[-\frac{1}{2\sigma^2}(x-x_0)^2]$ with $\sigma = 0.5$, $x_0 = (0.6, 0.25)$. The following figure shows the squared error $\|u_{\text{exact}}-u\|^2_{L_2}$ as a function of the polynomial degree $p$:

![Relative error](figures/error_reference_triangle.png)

Results are shown both for single precision and double precision arithmetic. We would expect that the error decreases for higher values of $p$ since the solution can be approximated better by higher degree polynomials. Although initially this is indeed the case, it appears that the error can not be reduced below a certain value and it in fact increases for larger values of $p$. To understand this behaviour, we need to discuss how (real) numbers are represented on a computer.

## Floating point numbers
A general **floating point number system** $\mathbb{F}$ is specified by four integer numbers:
* a base $1<\beta\in\mathbb{N}$
* a precision $0<p\in\mathbb{N}$
* a range of exponents defined by $L,U\in\mathbb{Z}$ with $L<0\le U$

The set $\mathbb{F}$ consists of all numbers $x$ of the form

$$
x = \pm \underbrace{\left(d_0 + d_1\beta^{-1} + d_2\beta^{-1} + \dots+d_{p-1}\beta^{1-p}\right)}_{\text{mantissa}}\cdot\beta^E\qquad(\dagger)
$$

where the coefficients $d_i\in \{0,1,2,\dots,\beta-1\}$ and the **exponent** $E$ with $L\le E\le U$ are natural numbers. The expression in brackets is called the **mantissa**. Note that although $\beta,p,L,U$ as well as $E,d_i$ are integers, they represent real numbers through $(\dagger)$.
  
The floating point number system $\mathbb{F}$ is called *normalised* if $d_0>0$; this makes each number in $\mathbb{F}$ unique.

#### Example
The number $234.7$ is
$$
\left(2+3\cdot 10^{-1}+4\cdot 10^{-2} +7\cdot 10^{-3}\right)\cdot 10^2
$$

in precision 4, base 10 arithmetic. It cannot be represented exactly in precision 3, base 10 arithmetic (and would have to be approximated as $\left(2+3\cdot 10^{-1}+5\cdot 10^{-2}\right)\cdot 10^2 = 235$ in this case).

The smallest positive normalised number of the from $(\dagger)$ is obtained by setting $d_0=1$, $d_i=0$ for $i>0$ and $E=L$. This results in $1\cdot \beta^L$ which is also is called the **underflow threshold**.

The largest positive normalised number in $\mathbb{F}$ is obtained by setting $d_i=\beta-1$, for $i\le 0$ and $E=U$. This results in

$$
\begin{aligned}
(\beta-1)\left(1+\beta^{-1}+\beta^{-2}+\dots+\beta^{1-p}\right)\cdot \beta^U &=
(\beta-1) \beta^U \sum_{j=0}^{p-1} \beta^{-p}\\
&= (1-\beta^{-p})\beta^{U+1},
\end{aligned}
$$

which is also called the **overflow threshold**.

Obviously, the floating point number system $\mathbb{F}$ is not closed under standard arithmetic operations: for example, $x,y\in\mathbb{F}$ does not necessarily imply that $x+y\in\mathbb{F}$. If a computation with two numbers $x,y\in\mathbb{F}$ results in a number $z\not\in\mathbb{F}$ we need to somehow represent $z$ by some nearby element $\tilde{z} := \mathcal{R}_{\mathbb{F}}(z)\in\mathbb{F}$ such that $|z-\widetilde{z}|$ is small. A common choice is to employ some sort of rounding.

### IEEE 754 (Normalised) Arithmetic
The most commonly used floating point systems on modern computers are single- and double-precision, which are implemented according to the [IEEE 754 standard](https://standards.ieee.org/ieee/754/6210/). In both cases, $\beta=2$ and it is implicitly assumed that $d_0=1$, so this number does not have to be stored.

#### Single precision
`np.float32`: One binary digit (=bit) is used to store the sign, 8 for the exponent and 23 for the mantissa $\Rightarrow$ 32 bits (4 bytes) in total.

![bits single precision](figures/bits_single_precision.svg)

* $p=24$
* $L=-126$, $U=127$
* Underfloat threshold = $2^{-126} \approx 10^{-38}$
* Overflow threshold = $2^{128} \approx 10^{38}$


#### Double precision
`np.float64`: One binary digit is used to store the sign, 11 for the exponent and 52 for the mantissa $\Rightarrow$ 64 bits (8 bytes) in total.

![bits double precision](figures/bits_double_precision.svg)

* $p=53$
* $L=-1022$, $U=1023$.
* Underfloat threshold = $2^{-1022} \approx 10^{-308}$
* Overflow threshold = $2^{1024} \approx 10^{308}$

### Representation of special values
Since $d_0=1$ it appears that we can not store the number zero. To represent this number and some other special cases, several dedicated bit patterns are reserved:

* The number zero is stored as `s000 0000 0000 0000 0000 0000 0000 0000` where $s$ is the sign bit. Note that there is both $+0$ ($s=1$) and $-0$ ($s=0$).
* `NaN` ("not a number") is stored as `s111 1111 1xxx xxxx xxxx xxxx xxxx xxxx` where the sequence denotes with $x$ stands for any non-zero number and the sign $s$ is usually ignored. The result of mathematically invalid operations (such as taking the square root of a negative number) is stored as a `NaN`
* Infinity ($\pm\infty$) is stored as `s111 1111 1000 0000 0000 0000 0000 0000` where again $s$ denotes the sign. This result can arise from division by zero.

### Machine epsilon
From $(\dagger)$ it can be seen that gaps between numbers in $\mathbb{F}$ increase for larger numbers. For each exponent $E$ the interval $[2^E,2^{E+1}]$ is discretised into $2^{p-1}$ equal pieces of size $2^{1-p}\cdot 2^E$, as shown in the following figure:

![floating point number spacing](figures/floating_point_spacing.svg)

Setting $E=0$, we see that the size of gap of numbers in $\mathbb{F}$ around $1$ is

$$
2^{1-p} = \begin{cases}
2^{-23} \approx 10^{-7} & \text{in single precision}\\
2^{-52} \approx 2\cdot 10^{-16} & \text{in double precision}
\end{cases}
$$

This quantity is also known as the **machine** epsilon $\varepsilon_{\text{mach}}$. It is the smallest positive number in $\mathbb{x}$ that can be added to $1$ such that (after rounding to $\mathbb{F}$) the result is different from $1$:

$$
\varepsilon_{\text{mach}} := \min_{x\in\mathbb{F},x>0}\{x: \mathcal{R}_{\mathbb{F}}(1+x)\neq 1\}
$$

Put differently, the machine epsilon is the relative size of rounding errors or the relative size of floating point operations. If $z$ is the result of some arithmetic operation involving numbers from $\mathbb{F}$ then

$$
\frac{|\mathcal{R}_{\mathbb{F}}(z)-z|}{|z|} \sim \varepsilon_{\text{mach}}.
$$

## Rounding errors
As the following examples show, rounding errors can have serious consequences

### Example 1 (harmless)
Consider the two numbers $x=4.7\cdot 10^{-16}$ and $x=2.9\cdot 10^{-16}$. Both can be represented exactly as floating point numbers. The same is true for their difference $z=x-y=1.8\cdot 10^{-16}$, i.e. $\widetilde{z}=\mathcal{G}_{\mathbb{F}}(z)=z$ and as a consequence the rounding error is zero of this operation is zero:

$$
\frac{\mathcal{R}_{\mathbb{F}}(z)-z}{z} = \frac{\mathcal{R}_{\mathbb{F}}(x)-\mathcal{R}_{\mathbb{F}}(y)-z}{z} = \frac{(4.7 -2.9 - 1.8)\cdot10^{-16}}{1.8\cdot 10^{-16}} = 0.
$$
In general, adding or subtracting numbers leads to small relative errors provided both the numbers and the result of the computation are of comparable size. As the following example shows, the final point is crucial.

### Example 2 (subtracting two numbers that are very close)
Now assume that we compute the difference of the two numbers by first adding $x$ and $y$ to one and then subtract the resulting numbers:
$x'=1+x$, $y'=1+y$, $z'=x'-y'$. Although in exact arithmetic $z'$ will be identical to $x-y$, this is not true in floating point arithmetic. First observe that $x'$ will be rounded to $\mathcal{R}_{\mathbb{F}}(x')=1.00000000000000044$ and $y'$ will be rounded to $\mathcal{R}_{\mathbb{F}}(x')=1.00000000000000022$. The relative error of $z'$ is:

$$
\frac{\mathcal{R}_{\mathbb{F}}(z')-z}{z} = \frac{\mathcal{R}_{\mathbb{F}}(x')-\mathcal{R}_{\mathbb{F}}(y')-z}{z} = \frac{(4.4 - 2.2 - 1.8)\cdot10^{-16}}{1.8\cdot 10^{-16}} = \frac{0.4}{1.8} \approx 23\%
$$

### Example 3 (adding two numbers of very different size)
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
and therefore $x_0=0$. Altogether we find $x_0=0$, $x_1=1$. This is very different from the exact solution $x_0=-\frac{1}{1+10^{-20}}$, $x_1=+\frac{1}{1+10^{-20}}$! Hence, although the rounding we performed in the numerical solution procedure appears to be innocent, we get a completely wrong solution. In this case, this problem can be fixed by using a slightly different solution procedure: subtract the first equation from the second equation to obtain $(1+10^{-20}) x_0 = -1$, which can be safely approximated by $x_0=-1$. Then use this in the second equation $x_0+x_1=0$ to conclude $x_1=+1$. Now this is very close to the exact solution of the perturbed system. It turns out that solving the perturbed linear system with [`np.linalg.solve()`](https://numpy.org/doc/2.2/reference/generated/numpy.linalg.solve.html) will also give a good approximation. Unfortunately, this is not always the case, as the following example demonstrates.

### Example 4 (ill-conditioned matrix)
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
In general, when solving linear systems with $n$ equations we are interested in quantifying the error on the solution.

In exact arithmetic the solution $\boldsymbol{u}\in\mathbb{R}^n$ satisfies

$$
A\boldsymbol{u} = \boldsymbol{b}
$$

where $A$ is a $n\times n$ matrix. However, due to rounding errors, the solution $\boldsymbol{u}'=\boldsymbol{u}+\delta\boldsymbol{u}$ that we actually compute corresponds to the system

$$
(A+\delta A)(\boldsymbol{u}+\delta\boldsymbol{u}) = \boldsymbol{b} + \delta\boldsymbol{b}
$$

If we set $\delta\boldsymbol{b}=0$ for the moment, we can derive a bound on $\delta\boldsymbol{u}$ for the case where Gaussian elimination is used to solve the linear system. For this we use the following norm on vectors and matrices:

$$
\begin{aligned}
\|\boldsymbol{w}\|_\infty &:= \max_{i=0,\dots,n-1} |w_i|,\\
\|A\|_\infty &:= \max_{\boldsymbol{w}\in\mathbb{R}^n,\boldsymbol{w}\neq \boldsymbol{0}} \frac{\|A\boldsymbol{w}\|_\infty}{\|\boldsymbol{w}\|_\infty} = \max_{0\le i < n} \sum_{j=0}^{n-1} |A_{ij}|
\end{aligned}
$$

Using the definition, it is easy to see that $\|A\boldsymbol{w}\|_\infty\le \|A\|_\infty\|\boldsymbol{w}\|_\infty$ and $\|AB\|_\infty\le\|A\|_\infty\|B\|_\infty$.

The condition number of a matrix $A$ is given by

$$
\text{cond}(A) := \|A^{-1}\|_\infty\|A\|_\infty.
$$

(For a real-valued symmetric positive matrix this is identical to the ratio of the largest and smallest eigenvalue).

We now have the following

### Theorem 1
If the $n\times n$ matrix $A$ is non-singular and $\delta A$ is sufficiently small, namely

$$
\|\delta A\|_{\infty} \|A^{-1}\|_\infty \le \frac{1}{2}
$$

then $A+\delta A$ is non-singular and

$$
\frac{\|\delta\boldsymbol{u}\|_\infty}{\|\boldsymbol{u}\|_\infty} \le 2\text{cond}(A) \frac{\|\delta A\|_\infty}{\|A\|_\infty}.
$$

Furthermore, it can be shown that if Gaussian elimination is used to solve the linear system (and this is the method used by [numpy.linalg.solve](https://numpy.org/doc/2.2/reference/generated/numpy.linalg.solve.html)), then the effect of roundoff errors is

$$
\frac{\|\delta A\|_\infty}{\|A\|_\infty} \le n \varepsilon_{\text{mach}} G(A),
$$

where $G(A)$ is a "well-behaved" number that depends on the matrix $A$. Although it is possible to construct pathological examples for which $G(A)$ is large, it is reasonable to assume that for the matrices that we consider here $G(A)$ is small and only depends weakly on $A$. Putting everything together, we find that

### Theorem 2
Under the conditions of Theorem 1 and if Gaussian elimination is used to solve the linear system, the error $\delta \boldsymbol{u}$ can be bounded

$$
\frac{\|\delta\boldsymbol{u}\|_\infty}{\|\boldsymbol{u}\|_\infty} \le 2G(A)\cdot n\cdot \text{cond}(A)\cdot \varepsilon_{\text{mach}}.
$$

### Summary
Although this is only an upper bound which does not necessarily have to be tight, this result implies that the (relative) error $\|\delta\boldsymbol{u}\|_\infty/\|\boldsymbol{u}\|_\infty$

* is proportional to the machine epsilon $\varepsilon_{\text{mach}}$,
* increases with problem size $n$,
* and grows with the condition number $\text{cond}(A)$ of the matrix.

### Estimating the error
In general it is not possible to compute the error $\delta\boldsymbol{u}$. However, since $\boldsymbol{b}-A\boldsymbol{u}=\boldsymbol{0}$, it is natural to consider the residual $\boldsymbol{r}:=\boldsymbol{b}-A\boldsymbol{u}'$ which measures to which degree the numerical solution $\boldsymbol{u}'$ fails to satisfy the linear system. Unfortunately, if $\|\boldsymbol{r}\|_\infty$ is small this does not necessarily imply the smallness of $\|\boldsymbol{u}\|_\infty$. To see this, first observe that $\boldsymbol{r}=A\delta\boldsymbol{u}$. Then note that since
$\|\boldsymbol{b}\|_\infty \le \|A\|_\infty \|\boldsymbol{u}\|_\infty$ and $\|\delta\boldsymbol{u}\|_\infty \le \|A^{-1}\|_\infty \|\boldsymbol{r}\|_\infty$

$$
\frac{\|\delta \boldsymbol{u}\|_\infty}{\|\boldsymbol{u}\|_\infty} 
\le \|A^{-1}\|_\infty\| \boldsymbol{r}\|_\infty \frac{\|A\|_\infty}{\|\boldsymbol{b}\|_\infty}\\
= \text{cond}(A)\frac{\|\boldsymbol{r}\|_\infty}{\|\boldsymbol{b}\|_\infty}
$$

Hence, the smallness of $\|\boldsymbol{r}\|_\infty/\|\boldsymbol{b}\|_\infty$ only implies the smallness of the relative error if the condition number of $A$ is small.

### Proof of Theorem 1 (not examinable)
Observe first that if $X$ is any $n\times n$ real-valued matrix with $\|X\|_\infty<1$ then $\|X^n\|_\infty\le \|X\|_\infty^n\rightarrow 0$ as $n\rightarrow\infty$. Thus

$$
(I-X)(1+X+X^2+\dots+X^n) = 1-X^{n+1} \rightarrow I\quad\text{as $n\rightarrow\infty$}.
$$

This implies that

$$
(I-X)^{-1}=\sum_{j=0}^{\infty} X^j
$$

and

$$
\|(I-X)^{-1}\|_\infty \le \sum_{j=0}^{\infty} \|X^j\|_\infty \le \sum_{j=0}^{\infty} \|X\|_\infty^j = (1-\|X\|_\infty)^{-1} \qquad(\star)
$$

Now write

$$
A + \delta A = (I+\delta A\;A^{-1})A
$$

and set $X=-\delta A\;A$. By the assumption $\|X\|_\infty\le \frac{1}{2}<1$ and we can apply $(\star)$ to show that $I+\delta A\;A^{-1}$ is non-singular and that

$$
\|(I+\delta A\;A)^{-1}\|_\infty\le (1-\|\delta A\;A^{-1}\|_\infty)^{-1} \le 2.
$$

As a consequence, $A+\delta A$ is non-singular and

$$
\|(A+\delta A)^{-1}\|_\infty = \|A^{-1}(I+\delta A\;A^{-1})^{-1}\|_\infty \le 2\|A^{-1}\|_\infty
$$

Finally, subtract the two equations $(A+\delta A)(\boldsymbol{u}+\delta\boldsymbol{u})=\boldsymbol{b}$ and $A\boldsymbol{u}=\boldsymbol{b}$ to obtain $(A+\delta A)\delta\boldsymbol{u} = -\delta A \boldsymbol{u}$. After multiplication by the inverse of $A+\delta A$ this becomes

$$
\delta\boldsymbol{u} = (A+\delta A)^{-1}\delta A\boldsymbol{u}.
$$

Taking the norm leads to

$$
\begin{aligned}
\|\delta \boldsymbol{u}\|_\infty &\le \|(A+\delta A)^{-1}\|_\infty \|\delta A\|_\infty \|\boldsymbol{u}\|_\infty\\
&\le 2 \underbrace{\|A^{-1}\|_\infty\|A\|_\infty}_{=\text{cond}(A)} \frac{\|\delta A\|_\infty}{\|A\|_\infty}\|\boldsymbol{u}\|_\infty
\end{aligned}
$$

To finish the proof, divide both sides of this inequality by $\|\boldsymbol{u}\|_\infty$.