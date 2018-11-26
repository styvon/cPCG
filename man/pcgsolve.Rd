\name{pcgsolve}
\alias{pcgsolve}
\alias{preconditioner}

\title{
Preconditioned conjugate gradient method
}
\description{
Preconditioned conjugate gradient method for solving system of linear equations Ax = b, where A is symmetric and positive definite, b is a column vector.  
}
\usage{
pcgsolve(A, b, preconditioner = "Jacobi", float tol = 1e-6, int maxIter = 1000)
}

\arguments{
  \item{A}{
  matrix, symmetric and positive definite.
  }
  \item{b}{
  vector, with same dimension as number of rows of A.
  }
  \item{preconditioner}{
  string, method for preconditioning: \code{"Jacobi"} (default), \code{"SSOR"}, or \code{"ICC"}.
  }
  \item{tol}{
  numeric, threshold for convergence, default is \code{1e-6}.
  }
  \item{maxIter}{
  numeric, maximum iteration, default is \code{1000}.
  }
}
\details{
When the condition number for \eqn{A} is large, the conjugate gradient (CG) method may fail to converge in a reasonable number of iterations. The Preconditioned Conjugate Gradient (PCG) Method applies a precondition matrix \eqn{C} and approaches the problem by solving: \deqn{{C}^{-1} A x = {C}^{-1} b} where the symmetric and positive-definite matrix \eqn{C} approximates \eqn{A} and \eqn{{C}^{-1} A } improves the condition number of \eqn{A}.  

Common choices for the preconditioner include: Jacobi preconditioning, symmetric successive over-relaxation (SSOR), and incomplete Cholesky factorization [2]. 

}

\section{Preconditioners }{
  \code{Jacobi}: The Jacobi preconditioner is the diagonal of the matrix A, with an assumption that all diagonal elements are non-zero.  
  
  \code{SSOR}: The symmetric successive over-relaxation preconditioner, implemented as \eqn{M = (D+L) D^{-1} (D+L)^T}. [1]  
  
  \code{ICC}: The incomplete Cholesky factorization preconditioner. [2]
}

\value{
Returns a vector representing solution x.
}
\references{
[1] David Young. “Iterative methods for solving partial difference equations of elliptic type”. In: Transactions of the American Mathematical Society 76.1 (1954), pp. 92–111.  
 
[2] David S Kershaw. “The incomplete Cholesky—conjugate gradient method for the iter- ative solution of systems of linear equations”. In: Journal of computational physics 26.1 (1978), pp. 43–65.  
}

\section{Warning }{
Users need to check that input matrix A is symmetric and positive definite before applying the function.  
}


\seealso{
\code{\link{cgsolve}}
}
\examples{
\dontrun{
test_A <- matrix(c(4,1,1,3), ncol = 2)
test_b <- matrix(1:2, ncol = 1)
pcgsolve(test_A, test_b, "ICC")
}
}

\keyword{ methods }
\keyword{ optimize }
