----

*&#169; Eike Mueller, University of Bath 2025. These notes are copyright of Eike Mueller, University of Bath. They are provided exclusively for educational purposes at the University and are to be downloaded or copied for your private study only. Further distribution, e.g. by upload to external repositories, is prohibited. html generated with [pandoc](https://pandoc.org/) using [easy-pandoc-templates](https://github.com/ryangrose/easy-pandoc-templates) under the [GPL-3.0.1 license](https://github.com/ryangrose/easy-pandoc-templates?tab=GPL-3.0-1-ov-file#readme)*

----
## Implementation

The implementation is straightforward: we construct a class `ThreePointQuadratureReferenceTriangle`, which is derived from the abstract base class `Quadrature` in [fem/quadrature.py](https://github.com/eikehmueller/finiteelements/blob/main/src/fem/quadrature.py). In the initialisation method, we set the nodes and associated weights:

```Python
self._nodes = np.asarray([[1 / 6, 1 / 6], [2 / 3, 1 / 6], [1 / 6, 2 / 3]])
self._weights = np.asarray([1 / 6, 1 / 6, 1 / 6])
```

The properties `nodes` and `weights` simply return these arrays:

```Python
@property
def nodes(self):
    return self._nodes

@property
def weights(self):
    return self._weights
```

The degree of precision is 2, and hence we return this in the `degree_of_precision()` method:

```Python
@property
def degree_of_precision(self):
    return 2
```

In the tests, we verify that the quadrature rule integrates monomials $x_0^{s_0}x_1^{s_1}$ of degree up to degree 2 exactly. `pytest.mark.parametrize()` is used to construct pairs of tuples $(s_0,s_1)$ and integrals $\int_{\widehat{K}} x_0^{s_0}x_1^{s_1}\;dx=\frac{s_0!s_1!}{(s_0+s_1+2)!}.$

The code can be found in [threepointquadrature.py](threepointquadrature.py) and the tests in [test_threepointquadrature.py](test_threepointquadrature.py).

## Theory
Using the given values of $w_q$ and $\xi^{(q)}$ we have that

$$
\begin{aligned}
S(s_0,s_1) &:= \sum_{q=0}^{2} w_q (\xi_0^{(q)})^{s_0}(\xi_1^{(q)})^{s_1}\\
&= \frac{1}{6}\left(\left(\frac{1}{6}\right)^{s_0}\left(\frac{1}{6}\right)^{s_0}+\left(\frac{2}{3}\right)^{s_0}\left(\frac{1}{6}\right)^{s_0}+\left(\frac{1}{6}\right)^{s_0}\left(\frac{2}{3}\right)^{s_0}\right)\\
&= \left(\frac{1}{6}\right)^{s_0+s_1+1}\left(1+4^{s_0}+4^{s_1}\right)
\end{aligned}
$$

Because of the symmetry $s_0\leftrightarrow s_1$ it is sufficient to consider four cases:
* $s_0=0$, $s_1=0$: $S(0,0) = \left(\frac{1}{6}\right)^{1}\left(1+4^{0}+4^{0}\right)=\frac{1}{2}=\frac{1}{2!} =\frac{0!0!}{(0+0+2)!}\\[1ex]$
* $s_0=1$, $s_1=0$: $S(1,0) = \left(\frac{1}{6}\right)^{2}\left(1+4^{1}+4^{0}\right)=\frac{1}{6}=\frac{1}{3!}=\frac{1!0!}{(1+0+2)!}\\[1ex]$
* $s_0=1$, $s_1=1$: $S(1,1) = \left(\frac{1}{6}\right)^{3}\left(1+4^{1}+4^{1}\right)=\frac{1}{24}=\frac{1}{4!}=\frac{1!1!}{(1+1+2)!}\\[1ex]$
* $s_0=2$, $s_1=0$: $S(2,0) = \left(\frac{1}{6}\right)^{3}\left(1+4^{2}+4^{0}\right)=\frac{1}{12}=\frac{2!}{4!}=\frac{2!0!}{(2+0+2)!}\\[1ex]$

Hence, the degree of precision of the quadrature rule is at least 2.

Polynomials of higher degree are not integrated exactly. Consider for example $s_0=3$, $s_1=0$. Then

$$
S(3,0) = \left(\frac{1}{6}\right)^4 \left(1+4^3+1\right) = \frac{66}{1296} = \frac{11}{216} \neq \frac{1}{20} = \frac{3!0!}{(3+0+2)!}
$$
