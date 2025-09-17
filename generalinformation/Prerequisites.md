<div align="center">
  <p style="font-size:32px;">MA32070 Prerequisites and background reading</p>
</div>

In the following we list some prior knowledge that is expected for this course together with some links to further information. Most of this should be familiar to you, and you will become more confident with writing Python code for solving problems in Scientific Computing throughout this course. However, you might benefit from brushing up on you Python coding skills and mathematical concepts that you have seen in other courses at the beginning of the semester.

----

*&#169; Eike Mueller, University of Bath 2025. These notes are copyright of Eike Mueller, University of Bath. They are provided exclusively for educational purposes at the University and are to be downloaded or copied for your private study only. Further distribution, e.g. by upload to external repositories, is prohibited. html generated with [pandoc](https://pandoc.org/) using [easy-pandoc-templates](https://github.com/ryangrose/easy-pandoc-templates) under the [GPL-3.0.1 license](https://github.com/ryangrose/easy-pandoc-templates?tab=GPL-3.0-1-ov-file#readme)*

----

# Mathematical background

## Numerical analysis
We will use several basic concepts that you have come across in your Numerical Analsys lecture, in particular:

* Interpolation of functions with [Lagrange polynomials](https://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html)
* Numerical quadrature with the [Gauss-Legendre method](https://mathworld.wolfram.com/Legendre-GaussQuadrature.html)

## The Finite Element method
You should be familiar with the fundamental ideas behind the finite element method, as taught in *"MA30366: Numerical solution of elliptic partial differential equations"*. However, the focus of the present unit is on *implementing* the finite element method, not on proofs and theory. We will recapitulate the key concepts at the beginning of this course, if anything looks unfamiliar please refer back to the lecture notes of MA30366 and/or consult the references below. 

## Linear Algebra
I expect you to be familiar with fundamental concepts of linear algebra: for example, you should know how to manipulate matrices and vectors, know what eigenvalues and eigenvectors are and understand under which conditions a matrix is invertible.

### Tensors
In the following we will often work with tensors, which are $d$-dimensional arrays. A tensor $T$ of rank  $d$ is an object which can be indexed with $d$ integers $i_0,i_1,\dots,i_{d-1}$ where $0\le i_k < s_k$ for $k=0,1,2,\dots,d-1$, i.e. we can write $T_{i_0,i_1,\dots,i_{d-1}}\in \mathbb{R}$ for the tensor elements. The list $[s_0,s_1,\dots,s_{d-1}]$ is the **shape** of the tensor. We are already familiar with two special cases:

* rank 0 tensors are scalars, i.e. real numbers $s\in\mathbb{R}$. In this case the shape is the empty list $[]$.
* rank 1 tensors are vectors $\boldsymbol{v}\in \mathbb{R}^n$ with elements $v_i$ for $0\le i< n$. The shape of a vector is the single number $[n]$, i.e. the dimension of the vector.
* rank 2 tensors are $n\times m$ matrices $A\in \mathbb{R}^{n\times m}$ with elements $A_{ij}$ for $0\le i<n$ and $0\le j <m$. The shape is the tuple $[n,m]$, where $n$ is the number of rows and $m$ is the number of columns of the matrix.

#### Adding tensors
Tensors $T$, $T'$ of the same shape can be scaled and added: if $\alpha,\beta\in \mathbb{R}$ are real numbers, then $S=\alpha T+\beta T'$ is a new tensor of the same shape. The entries of $S$ are given by

$$
S_{i_0,i_1,\dots,i_{d-1}} = \alpha T_{i_0,i_1,\dots,i_{d-1}}+\beta T'_{i_0,i_1,\dots,i_{d-1}}
\qquad \text{for all $i_0,i_1,\dots,i_{d-1}$ with $0\le i_k < s_k$ for $k=0,1,2,\dots,d-1$}.
$$

#### Multiplying tensors
Tensors can be multiplied in different ways by contracting indices. For example, we can multiply the rank 3 tensor $T$ and the rank 4 tensor $T'$ to obtain a rank 5 tensor $R$ as follows:

$$
R_{ijm\ell n} = \sum_{k} T_{ijk} T'_{mk\ell n}
$$

Alternatively, we could also compute a rank 1 tensor $S$ from $T$, $T'$ and the rank 2 tensor $T''$:

$$
S_{\ell} = \sum_{ijk} T_{ijk} T'_{kjni} T''_{n\ell}
$$

Note that the order of the indices matters! The matrix-vector product $\boldsymbol{v}=A\boldsymbol{w}$ is a special case:

$$
v_i = \sum_j A_{ij}w_j
$$

Contractions can also result in a scalar. For example

$$
\boldsymbol{v}\cdot \boldsymbol{w} = \sum_i v_i w_i
$$

is the dot-product of two vectors $\boldsymbol{v}$ and $\boldsymbol{w}$, whereas

$$
\operatorname{trace}(A) = \sum_{i} A_{ii}
$$

is the trace of a matrix $A$.

# Python
To complete this course successfully you need to have a solid command of the Python language, which you should be familiar with from your first year *"Programming and Mathematics"* unit. Below you can find a (non-exhaustive) list of fundamental concepts that we will use frequently. If any of these ideas are new to you, follow the links or consult the [Python tutorial](https://docs.python.org/3/tutorial/index.html), [Python language reference](https://docs.python.org/3/reference/index.html) and the [numpy API reference](https://numpy.org/doc/stable/reference/index.html), in particular the section on [linear algebra](https://numpy.org/doc/stable/reference/routines.linalg.html). If you have never seen some of the concepts don't worry, since they will become more familiar over the course of this semester. A crucial skill of a successful scientific programmer is the ability to find and digest documentation: if you come across a concept that you are not familiar with you should actively try to understand it. For this, you might want to implement some simplified examples to see what is going on or use google to find additional documentation.

## Fundamental Python concepts

### Control flow and code structure
* Definining [functions](https://docs.python.org/3/tutorial/controlflow.html#defining-functions) and [passing arguments](https://docs.python.org/3/tutorial/controlflow.html#more-on-defining-functions)
* [for-loops](https://docs.python.org/3/tutorial/controlflow.html#for-statements), and the [range()](https://docs.python.org/3/tutorial/controlflow.html#the-range-function) function
* Control-flow with [if-statements](https://docs.python.org/3/tutorial/controlflow.html#if-statements)
* Organising functionality in [modules](https://docs.python.org/3/tutorial/modules.html)

### Data types
* integer and floating point [numbers](https://docs.python.org/3/tutorial/introduction.html#numbers)
* [strings](https://docs.python.org/3/tutorial/introduction.html#text) for storing characters and text
* [lists](https://docs.python.org/3/tutorial/introduction.html#lists) and [sets](https://docs.python.org/3/tutorial/datastructures.html#sets) to store collections of items
* [dictionaries](https://docs.python.org/3/tutorial/datastructures.html#dictionaries) to efficiently store (key $\rightarrow$ value) pairs.

### Object Oriented programming
* Defining [classes](https://docs.python.org/3/tutorial/classes.html)
* Class [methods](https://docs.python.org/3/tutorial/classes.html#method-objects) and [properties](https://docs.python.org/3/tutorial/classes.html#instance-objects)
* Reusing functionality with [inheritance](https://docs.python.org/3/tutorial/classes.html#inheritance)
* [Abstract classes](https://docs.python.org/3/library/abc.html)

### More advanced Python concepts
* Anonymous functions with [lambda expressions](https://docs.python.org/3/tutorial/controlflow.html#lambda-expressions)
* Manipulating lists with [list comprehensions](https://docs.python.org/3/tutorial/datastructures.html#list-comprehensions)
* Fancy [output formatting](https://docs.python.org/3/tutorial/inputoutput.html#fancier-output-formatting)
* [Decorators](https://peps.python.org/pep-0318/) for modifying methods

### Numpy
* Constructing [arrays](https://numpy.org/doc/stable/user/basics.creation.html)
* [Indexing and slicing](https://numpy.org/doc/stable/user/basics.indexing.html) arrays
* Matrix- and vector manipulation such as [matrix transposition](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.transpose.html) and [matrix multiplication](https://numpy.org/devdocs/reference/generated/numpy.matmul.html)
* Other basic [numerical linear algebra](https://numpy.org/doc/stable/reference/routines.linalg.html#module-numpy.linalg) operations

#### Tensors in numpy
In numpy, tensors are represented by multidimensional arrays of type [`np.ndarray(...,dtype=float)`](https://numpy.org/devdocs/reference/generated/numpy.ndarray.html). For example, we can create a rank 3 tensor of shape $(2,3,4)$ with only zero entries like this:

```Python
T = np.zeros(shape=(2,3,4),dtype=float)
```

or construct the $2\times 3$ matrix $\begin{pmatrix}1.8 & 2.2 & 3.4 \\ 4.2 & 5.1 & 6.7\end{pmatrix}$ like this:

```Python
A = np.asarray([[1.8,2.2,3.4],[4.2,5.1,6.7]],dtype=float)
```

The shape is given by `T.shape` and `A.shape`.

Tensors of the same shape can be scaled and added like this:

```Python
S = alpha*T + beta*Tprime
```

In fact, Python allows the addition of tensors of different shapes according to a set of [broadcasting rules](https://numpy.org/doc/stable/user/basics.broadcasting.html).

Tensors can be multiplied with the powerful [`numpy.einsum()` function](https://numpy.org/doc/stable/reference/generated/numpy.einsum.html). For example, to compute $R_{ijm\ell n} = \sum_{k} T_{ijk} T'_{mk\ell n}$ use

```Python
R = np.einsum("ijk,mkln->ijmln",T,Tprime)
```

and to compute $S_{\ell} = \sum_{ijk} T_{ijk} T'_{kjni} T''_{n\ell}$ use

```Python
S = np.einsum("ijk,kjni,nl->l",T,Tprime,Tprimeprime)
```

The matrix-vector product $\boldsymbol{v}=A\boldsymbol{w}$ can be written as

```Python
v = np.einsum("ij,j->i",A,w)
```

or alternatively simply as

```Python
v = A @ w
```

The latter us usually preferred since it is easier to understand what the code does.

Look at [the documentation of `numpy.einsum()`](https://numpy.org/doc/stable/reference/generated/numpy.einsum.html) for further details.

## Documentation and style
Comment your code as much as possible (but be concise). In particular, use [docstrings](https://docs.python.org/3/tutorial/controlflow.html#documentation-strings) to describe the input, output and functionality of functions. Use a consistent style and naming convention (for example, start class names with capital letters and variables with lowercase letters). You might want to take some inspiration form the [PEP8 style guide](https://peps.python.org/pep-0008/), but there is no need to strictly conform to this. It is much more important that you use a consistent style. VSCode (see below) can help with enforcing this. [Section 4 of [Ham23]](https://object-oriented-python.github.io/4_style.html) contains a useful discussion of style.

# Other tools

## Visual Studio Code
[Visual Studio (VS) Code](https://code.visualstudio.com/) is an extremely powerful editor which is used by many developers today. It supports static code analysis via [linting](https://code.visualstudio.com/docs/python/linting) and automatic [formatting](https://code.visualstudio.com/docs/python/formatting) which can be incredibly helpful. Furthermore, it comes with a built-in [debugger](https://code.visualstudio.com/docs/debugtest/debugging). You are free to use a different editor for this course, but it is quite difficult to beat VSCode's functionality.

## Version control with git
The de-facto standard for version control nowadays is [git](https://git-scm.com/). This is a tool which allows the easy tracking of changes in source code and collaboration with other developers through hosting sites such as [github](https://github.com/). We won't have time to discuss git in this course and it is not compulsory for you to use version control when developing your code, but it is recommended to make your life easier. A good introduction to git is the Software Carpentry course on [Version Control with git](https://swcarpentry.github.io/git-novice/). Here we will only need the git command line tool download the provided finite library. There are several graphical user interfaces such as [Sourcetree](https://www.sourcetreeapp.com/) and [GitKraken](https://www.gitkraken.com/) which you might want to try out; VSCode also has built-in git support.

## Testing
Remember that **all untested code is wrong.** When writing any code it is important to simultaneously implement tests that verify its correctness. Fortunately, this is very straightforward with the [pytest](https://docs.pytest.org/en/stable/) package.

Let's look at a simply example. To test a function `addition(a,b)` which adds the two numbers $a$ and $b$ we simply implement a function

```Python
def test_addition():
    assert addition(3,4) == 7
```

in a script called `test_summation.py`. Running `pytest -v` in the directory which contains `test_summation.py` will automatically look for all methods starting with `test_` and check whether the corresponding assertions are satisfied or not:

```
$ pytest -v

============================= test session starts ==============================
platform darwin -- Python 3.13.3, pytest-8.3.5, pluggy-1.5.0 -- /opt/homebrew/opt/python@3.13/bin/python3.13
cachedir: .pytest_cache
rootdir: /Users/eikehmueller/Work/Bath/git_workspace/ma32070
collected 1 item

test_summation.py::test_addition PASSED                                  [100%]

============================== 1 passed in 0.01s ===============================
```

## Debugging
Having a good strategy for debugging code is equally important, and you should already know the basics from your first year programming course. [Section 8 of [Ham23]](https://object-oriented-python.github.io/8_debugging.html) contains a detailed discussion of this topic and a description of some popular debugging tools, including the built-in debugger of VSCode. While you will notice that your debugging skills naturally improve through practice, and you should use the techniques and tools that best work for you, I encourage you to read David Ham's discussion of [Debugging Strategy](https://object-oriented-python.github.io/8_debugging.html#debugging-strategy).

## References
### Finite element methods
The following references might be useful for further reading on finite element methods. You are not expected to study this material in detail, but they might be helpful as a reference and if you want to deepen your understanding further.

* [[CH24]](https://finite-element.github.io/) Lecture notes on *"Finite elements: Analysis and implementation"* by Colin Cotter and David Ham (Imperial College)
* [[Far21]](https://people.maths.ox.ac.uk/farrellp/femvideos/notes.pdf) Lecture notes on  *"Finite Element Methods for PDEs"* by Patrick Farrell (Oxford)
* [[Log11]](http://launchpad.net/fenics-book/trunk/final/+download/fenics-book-2011-10-27-final.pdf) *"Automated Solution of Differential Equations by the Finite Element Method"* (the FEniCS book) by Anders Logg et al.

### Python
The following two books contain additional information on Python and how to use it to implement mathematical software.

* [[Ram22]](https://www.oreilly.com/library/view/fluent-python-2nd/9781492056348/) *"Fluent Python"* by Luciano Ramalho, 2nd edition 2022 (available through the University library)
* [[Ham23]](https://object-oriented-python.github.io/) *"Object-oriented Programming in Python for Mathematicians"* by David Ham

[[Ram22]](https://www.oreilly.com/library/view/fluent-python-2nd/9781492056348/) provides a comprehensive overview of Python concepts, including advanced topics. It is a useful reference if you want to develop a deeper understanding of the inner workings of Python. [[Ham23]](https://object-oriented-python.github.io/) uses Python as a tool to implement mathematical objects. It contains many higher-level ideas such as the use of suitable abstractions for writing well-designed code.


## Getting help
You might want to review the [Appendix 1. Getting Help](https://object-oriented-python.github.io/a1_help.html) from [David Ham's excellent book](https://object-oriented-python.github.io/index.html) on how to ask good questions.