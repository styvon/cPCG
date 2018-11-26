\name{cPCG-package}
\alias{cPCG-package}
\alias{cPCG}
\docType{package}
\title{
\packageTitle{cPCG}
}
\description{
\packageDescription{cPCG}  


}
\details{

Functions in this package serve the purpose of solving for x in \eqn{Ax = b}, where A is a symmetric and positive definite matrix, b is a column vector. 

To improve scalability of conjugate gradient methods for larger matrices, the \code{Armadillo} templated C++ linear algebra library is used for the implementation. The package also provides flexibility to have user-specified preconditioner options to cater for different optimization needs.  

The DESCRIPTION file:
\packageDESCRIPTION{cPCG}
\packageIndices{cPCG}

}
\author{
\packageAuthor{cPCG}
}
\references{
[1] Reeves Fletcher and Colin M Reeves. “Function minimization by conjugate gradients”. In: The computer journal 7.2 (1964), pp. 149–154.  

[2] David S Kershaw. “The incomplete Cholesky—conjugate gradient method for the iter- ative solution of systems of linear equations”. In: Journal of computational physics 26.1 (1978), pp. 43–65.   

[3] Yousef Saad. Iterative methods for sparse linear systems. Vol. 82. siam, 2003.  

[4] David Young. “Iterative methods for solving partial difference equations of elliptic type”. In: Transactions of the American Mathematical Society 76.1 (1954), pp. 92–111.  
}
\keyword{ package }

\examples{
  # generate test data
  test_A <- matrix(c(4,1,1,3), ncol = 2)
  test_b <- matrix(1:2, ncol = 1)
  
  # conjugate gradient method solver
  cgsolve(test_A, test_b, 1e-6, 1000)
  
  # preconditioned conjugate gradient method solver, 
  # with incomplete Cholesky factorization as preconditioner
  pcgsolve(test_A, test_b, "ICC")
}
