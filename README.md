# cPCG: Efficient and Customized Preconditioned Conjugate Gradient Method

An R package to solve system of linear equations using (preconditioned) conjugate gradient algorithm, with improved efficiency using Armadillo templated C++ linear algebra library, and flexibility for userspecified preconditioning method.


# Installation options
Download `cPCG_1.0.tar.gz` file and build from command line:
```
R CMD INSTALL cPCG_1.0.tar.gz
```

Get current development version from github:

```R
# install.packages("devtools")
devtools::install_github("styvon/cPCG")
```

# Functions

## `cgsolve`
Conjugate gradient method for solving system of linear equations Ax = b, where A is symmetric and positive definite, b is a column vector.  

```
cgsolve(A, b, float tol = 1e-6, int maxIter = 1000)
```

## `pcgsolve`

When the condition number for A is large, the conjugate gradient (CG) method may fail to converge in a reasonable number of iterations. The Preconditioned Conjugate Gradient (PCG) Method applies a precondition matrix C and approaches the problem by solving C^{-1} A x = {C}^{-1} b where the symmetric and positive-definite matrix C approximates A and C^{-1}A  improves the condition number of A.  

```
pcgsolve(A, b, preconditioner = "Jacobi", float tol = 1e-6, int maxIter = 1000)
```

Common choices for the preconditioner include: Jacobi preconditioning, symmetric successive over-relaxation (SSOR), and incomplete Cholesky factorization. 

* `Jacobi`: The Jacobi preconditioner is the diagonal of the matrix A, with an assumption that all diagonal elements are non-zero.  
  
* `SSOR`: The symmetric successive over-relaxation preconditioner, implemented as M = (D+L) D^{-1} (D+L)^T.  
  
* `ICC`: The incomplete Cholesky factorization preconditioner. 


```
Roxygen: list(roclets = c("rd", "collate"))
```
