# Numerical quadrature
## Gauss-Legendre quadrature in one dimension
Numerical quadrature aims to approximate the integral of a function with a finite sum:

$$
\int_{-1}^{+1} f(z)\;dz \approx \sum_{i=0}^{n-1} w_i f(z^{(i)})
$$

A particular quadrature rule $\{(z^{(i)},w_i)\}_{i=0}^{n-1}$ is defined by the sets of points $z^{(i)}$ and corresponding weights $w_i$. Here we will consider Gauss-Legendre quadrature, for which the points are the roots of the [Legendre polynomial](https://mathworld.wolfram.com/LegendrePolynomial.html) $P_n(z)$ and the weights are given by $w_i = \frac{2}{(1-(z^{(i)})^2)(P_n'(z^{(i)}))^2}$. The points and weights can be constructed with [numpy.polynomial.legendre.leggauss](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.legendre.leggauss.html):

```
points, weights = numpy.polynomial.legendre.leggauss(n)
```
Crucially, the numerical quadrature is exact if the function to be integrated is a polynomial of degree $2n-1$:

$$
\int_{-1}^{+1} p(z)\;dz = \sum_{i=0}^{n-1} w_i p(z^{(i)})\qquad\text{for $p\in\mathcal{P}_{2n-1}$}
$$
## Two-dimensional quadrature for the reference triangle
To numerically integrate functions over the reference triangle $K$, first observe that $K$ is the image of the square $Q=[-1,+1]\times [-1,+1]$ under the Duffy transform $\tau$ which maps a point $\widetilde{x}=(\widetilde{x}_0,\widetilde{x}_1)\in Q$ to

$$
\begin{aligned}
\tau(\widetilde{x}) = \left(\frac{1}{2}(1+\widetilde{x}_0),\frac{1}{4}(1-\widetilde{x}_0)(1+\widetilde{x}_1)\right)^\top \in K
\end{aligned}
$$

### Integration over $\boldsymbol{Q}$
Since $Q$ is the product of two intervals, we can integrate functions over $Q$ by picking quadrature points and weights $\{(\widetilde{\zeta}^{(i)},\widetilde{w}_i)\}_{i=0}^{n(n+1)-1}$ with 

$$
\widetilde{\zeta}^{(i)} = \left(\widetilde{\zeta}^{(j)}_0,\widetilde{\zeta}^{(k)}_1\right)^\top,\quad \widetilde{w}_i = \widetilde{w}_{0,j}\cdot \widetilde{w}_{1,k} \qquad \text{where $i=nj+k$}.
$$

Here $\{(\widetilde{\zeta}^{(j)}_0,\widetilde{w}_{0,j})\}_{j=0}^{n-1}$ are $\{(\widetilde{\zeta}^{(k)}_0,\widetilde{w}_{0,k})\}_{k=0}^{n}$ with $n+1$ and $n$ points respectively (we need to integrate more accurately in the $0$-direction since and additional factor of $\widetilde{x}_0$ is introduced by the Duffy-transform). 

### Integration over $\boldsymbol{K}$
The quadrature rule $\{(\zeta^{(i)},w_i)\}_{i=0}^{n(n+1)-1}$ over $K$ is then obtained as

$$
\begin{aligned}
\zeta^{(i)} &= \tau(\widetilde{\zeta}^{(i)}) = \left(\frac{1}{2}(1+\widetilde{\zeta}^{(j)}_0),\frac{1}{4}(1-\widetilde{\zeta}^{(j)}_0)(1+\widetilde{\zeta}^{(k)}_1)\right)^\top,\\
w_i &= \widetilde{w}_i \left|\det\left(\frac{\partial \tau}{\partial \widetilde{x}}\right)\right| = \frac{1}{8}\widetilde{w}_{0,1}\widetilde{w}_{0,1}(1-\widetilde{\zeta}^{(j)}_0)
 \qquad \text{where $i=nj+k$.}
\end{aligned}
$$

The following figure shows the quadrature points on $Q$ and $K$ for $n=2$ points.

![Quadrature points on $Q$ and $K$](figures/quadrature.png)

## Implementation in Python