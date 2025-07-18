# Solvers and Preconditioners
We now consider methods for solving the linear system $Au=b$ with PETSc. Although PETSc also supports direct solvers such as Gaussian elimination, the design philosophy of the library assumes that all solvers are iterative methods. This means that, starting from a starting guess, $u^{(0)}$, the solution is updated iteratively as $u^{(k+1)}\gets u^{(k)}$ such that (hopefully) $u^{(k)}$ converges to the true solution $u$. Direct solvers can be considered as a special case in which this iteration converges in a single step.

## PETSc solver architecture
PETSc solvers are separated into two components:

1. an *iterative solver* or `KSP` object which describes how to perform the update $u^{(k+1)}\gets u^{(k)}$; we have already seen this in the previous lecture.
2. a *preconditioner* or `PC` object which accelerates the iteration.

To motivate this architecture and explain what a preconditioner is, let us consider a very simple iterative method. For this, we multiply the linear equation by a matrix $P^{-1}$ to obtain the equivalent system

$$
P^{-1} A u = P^{-1} b
$$

The matrix $P$, which we assume to be full rank and invertible, is the preconditioner. More precisely, in PETSc a `PC` object provides a way of inverting $P$, i.e. solving $Pz=r$ for a given vector $r$. In principle, we could choose any matrix such as $P=\mathbb{I}$ the identity matrix or the $P=D$ the diagonal of $A$. A good preconditioner has two properties:

1. It should be "close" to $A$, i.e. $P\approx A$ (in some sense)
2. Multiplication by $P^{-1}$ should be inexpensive, i.e. solving the linear system $Pz=r$ for a given vector $r$ should be cheap.

## Richardson iteration with Jacobi preconditioner
As an example, we now construct a very simple iterative procedure. Assume that we know the approximate solution $u^{(k)}$. If we knew the error $e^{(k)} := u-u^{(k)}$, it would be possible to compute the solution by simply adding $e^{(k)}$ to $u^{(k)}$. Unfortunately, this is not possible since knowledge of $e^{(k)}$ would require knowledge of the exact solution! Instead, let's try to construct an approximation $r^{(k)}$ to the true error, namely
$$
r^{(k)} = P^{-1} A e^{(k)}
$$
If $P\approx A$, this will be a good approximation. Interestingly we can compute $\widetilde{e}^{(k)}$ since
$$
\begin{aligned}
r^{(k)} &= P^{-1} A (u-u^{(k)})\\
&= P^{-1}(Au-Au^{(k)})\\
&= P^{-1}(b-Au^{(k)})
\end{aligned}
$$
The quantity $r^{(k)}$ is also known as the preconditioned residual, since $b-Au^{(k)}$ measures by how much the approximate solution violates the equation $Au=b$.

This leads to the preconditioned Richardson iteration:
$$
\begin{aligned}
u^{(k+1)} &= u^{(k)} + r^{(k)}\\
&= u^{(k)} + P^{-1}(b-Au^{(k)})
\end{aligned}
$$
It can be shown that
$$
e^{(k+1)} = \left(\mathbb{I} - P^{-1} A\right)e^{(k)}
$$
Hence, if the spectral radius of $\mathbb{I} - P^{-1} A$ is smaller than $1$, the iteration will converge to the true solution.

Note in particular that if we set $P=A$, the iteration converges in a single step.

One possible preconditioner is the diagonal $D$ of $A$, i.e.
$$
D_{ij} = \begin{cases}
A_{ii} & \text{if $i=j$}\\
0 & \text{otherwise}
\end{cases}
$$
This matrix is very simple to invert: to solve $Pz=r$ we can compute $z_i = r_i A_{ii}$. The choice $P=D$ is also know as the *Jacobi* preconditioner.

### PETSc options
To use a particular solver/preconditioner combination in PETSc, we need to specify solver options via the command line. For this, we need to add the following at the beginning of our Python script:
```python
import sys
import petsc4py

petsc4py.init(sys.argv)
```
This will ensure than any options passed to the script are parsed by PETSc. Then, when setting up the `KSP` object, we pass these options by calling 

```python
ksp.setFromOptions()
```
For example, to use the Richardson iteration with Jacobi preconditioner, we call our script like this:

```bash
python script.py -ksp_type richardson -pc_type jacobi
```
To monitor convergence, we can add the option `-ksp_monitor`, which will print out the norm $\|r^{(k)}\|$ of the (preconditioned) residual at each iteration. The iteration will stop once the residual norm has been reduced by a factor of at least $\epsilon$, i.e. $\|r^{(k)}\|/\|r^{(0)}|\|<\epsilon$. This tolerance can be controlled by `-ksp_rtol epsilon` (it is also possible to set an absolute convergence criterion $\|r^{(k)}\|<\epsilon_{\text{abs}}$ with `-ksp_atol`). Furthermore, we can tell PETSc to print information on the solver to some file `ksp_view.txt` with `-ksp_view :ksp_view.txt`. So, in summary to solve to a relative tolerance of $\epsilon=10^{-9}$ with the Jacobi-preconditioned Richardson iteration we would call:

```bash
python script.py -ksp_type richardson -pc_type jacobi -ksp_monitor -ksp_rtol 1.0E-9 -ksp_view :ksp_view.txt
```

The output looks like this:
```
  0 KSP Residual norm 1.838591537060e+00
  1 KSP Residual norm 8.624966321088e-01
  2 KSP Residual norm 3.904353664205e-02
  3 KSP Residual norm 1.230500385292e-02
  4 KSP Residual norm 1.919611985913e-03
  5 KSP Residual norm 6.049870199858e-04
  6 KSP Residual norm 9.437951818362e-05
  7 KSP Residual norm 2.974475251910e-05
  8 KSP Residual norm 4.640257259221e-06
  9 KSP Residual norm 1.462428569925e-06
 10 KSP Residual norm 2.281425870008e-07
 11 KSP Residual norm 7.190166811597e-08
 12 KSP Residual norm 1.121684357597e-08
 13 KSP Residual norm 3.535112890997e-09
 14 KSP Residual norm 5.514865252918e-10
```

and we can inspect the file `ksp_view.txt` to double check that the solver options have been set correctly:

At the beginning, it lists details on the `KSP Object`, which describes the iterative solver:

```
KSP Object: 1 MPI process
  type: richardson
    damping factor=1.
  maximum iterations=10000, initial guess is zero
  tolerances: relative=1e-09, absolute=1e-50, divergence=10000.
  left preconditioning
  using PRECONDITIONED norm type for convergence test
```

This is followed by details on the `PC Object`, which describes the preconditioner:
```
PC Object: 1 MPI process
  type: jacobi
    type DIAGONAL
  linear system matrix = precond matrix:
  Mat Object: 1 MPI process
    type: seqdense
    rows=5, cols=5
    total: nonzeros=25, allocated nonzeros=25
    total number of mallocs used during MatSetValues calls=0
```
#### Exercise
Instead of $P=D$, we could also use the lower triangular part of $A$ and set $P=D+L$ where
$$
L = \begin{cases}
A_{ij} & \text{if $i<j$} \\
0 & \text{otherwise}
\end{cases}
$$
Convince yourself that for a given vector $r$ the equation $(D+L)z=r$ can be solved row-by-row, i.e. by computing first $z_0 = r_0/A_{00}$, then computing $z_1 = (r_1 - A_{10}z_0)/A_{11}$, $z_2=(r_2 - A_{20}z_0 - A_{21}z_1)/A_{22}$ and so once. The corresponding preconditioner is also known as the successive overrelaxation (SOR) method. It can be chosen by setting `-pc_type sor``. Run the code with this preconditioner - how does the number of iterations change?

### Direct solvers
PETSc also supports direct solvers, which are implemented as preconditioners. For example, to use Gaussian elimination, we would set `-pc_type lu`. In this case, PETSc computes the factorisation $P=A=LU$, where $L$ and $U$ and lower- and upper-triangular matrices. Knowing $L$ and $U$ we can solve $Az=LUz=r$ by solving $Lz'=r$ and then $Uz=z'$. In this case, iterative solver will converge in a single iteration:
```
  0 KSP Residual norm 2.291531974978e+00
  1 KSP Residual norm 2.640425190731e-16
```
In this case, we can request that PETSc only applies the preconditioner, i.e. computes $u^{(1)} = P^{-1}b$ directly. Be careful with using `-ksp_type preonly`: if the preconditioner is not a direct solver, the iteration will simply stop after one iteration and return an incorrect result. For example, `-ksp_type preonly -pc_type richardson` will print out
```
  0 KSP Residual norm 1.405809375413e+01
  1 KSP Residual norm 6.204942588128e+00
```
and it is up to us to recognise that the computed solution does not solve $Au=b$.
