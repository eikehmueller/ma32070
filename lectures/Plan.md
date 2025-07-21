# Plan 

## 0. Background reading
* command line
* git and github (classroom?)
* Advanced Python concepts
  * Classes and inheritance
  * Abstract base classes
  * numpy
  * generators
  * properties
* VSCode
  * connecting to the cluster
  * formatting and linting

## 1. Motivation and background
* Introduction and housekeeping
  * Background reading
  * Structure of the course
  * Assessment
* Mathematical background
  * The Poisson-like problem
  * Weak formulation
  * Finite element subspaces
  * The resulting linear algebra problem
  * Solvers
  * Theory: existing of solution, convergence $\Rightarrow$ *``Numerical Solultion of Elliptic Problems''*

## 2. Finite elements
* Consider a domain consisting of a single triangular cell
  * Numbering of vertices and edges
* Ciarlet's definition of the finite element
  * Discuss basis of dual space and nodes
* Lagrange elements
  * Association with topological entities
  * Basis functions
  * Vandermonde matrix
  * Tabulating finite elements: values and derivatives
* Implementation in Python

## 3. Quadrature and solver on a single triangle
* Quadrature
  * 1d quadrature
  * extension to 2d
* Algorithms on a single element
  * Interpolating functions
  * Matrix assembly
* Python implementation

## 4. Interlude: rounding errors
* Numerical investigation of single-element solver
* Solution of a simple $2\times 2$ system
* Backward error analysis 

## 5. (Unstructured) meshes
* Topological entities and (co-) dimensions
* Encoding connectivity information
* Grid cells and reference elements
* Encoding geometry information: curved elements
* Implementation in Python

## 6. Function spaces
* Local-to-global maps
* Implementation in Python

## 7. Assembling and solving the global problem
* Assembling the matrix
* Assembling the right-hand side
* Application to test problem: manufactured solutions
  * Solution of the linear problem
  * Error and convergence
* General discussion of errors
  * Modelling error
  * Discretisation error
  * Rounding error / machine precision

## 8. Sparse matrices
* Compressed sparse row storage
* Introduction to PETSc and petsc4py
  * Key ideas
  * Fundamental data structures: Mat, Vec, KSP
* Assembling into PETSc matrices and vectors
* Implementation in Python

## 9. Solvers and preconditioners
* Iterative solvers
  * General idea (compare to direct solvers)
  * Krylov subspace methods
  * Preconditioning: Jacobi, multigrid, ILU
* Implementation in Python and petsc4py
* Numerical experiments
  * Comparison of number of iterations
  * Robustness

## 10. Firedrake
* General concepts, programming models, performance indicators
* Implementation of Poisson-like problem

## 11. Parallelisation
* Parallelisation strategies
* Domain decomposition
* Demo in Firedrake

## X. Advanced concepts
* Mixed problems
* DG methods
* Implementing Dirichlet boundary conditions
* Non-linear problems

