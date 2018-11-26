\name{cgsolve}
\alias{cgsolve}

\title{
Conjugate gradient method
}
\description{
Conjugate gradient method for solving system of linear equations Ax = b, where A is symmetric and positive definite, b is a column vector.  
}
\usage{
cgsolve(A, b, float tol = 1e-6, int maxIter = 1000)
}

\arguments{
  \item{A}{
  matrix, symmetric and positive definite.
  }
  \item{b}{
  vector, with same dimension as number of rows of A.
  }
  \item{tol}{
  numeric, threshold for convergence, default is \code{1e-6}.
  }
  \item{maxIter}{
  numeric, maximum iteration, default is \code{1000}.
  }
}
\details{
The idea of conjugate gradient method is to find a set of mutually conjugate directions for the unconstrained problem \deqn{ arg min_x f(x)} where \eqn{f(x) = 0.5 b^T A b - bx + z} and \eqn{z} is a constant. The problem is equivalent to solving \eqn{Ax = b}.  

This function implements an iterative procedure to reduce the number of matrix-vector multiplications [1]. The conjugate gradient method improves memory efficiency and computational complexity, especially when \eqn{A} is relatively sparse. 
}
\value{
Returns a vector representing solution x.
}
\references{
[1] Yousef Saad. Iterative methods for sparse linear systems. Vol. 82. siam, 2003.
}

\section{Warning }{
Users need to check that input matrix A is symmetric and positive definite before applying the function.  
}


\seealso{
\code{\link{pcgsolve}}
}
\examples{
\dontrun{
test_A <- matrix(c(4,1,1,3), ncol = 2)
test_b <- matrix(1:2, ncol = 1)
cgsolve(test_A, test_b, 1e-6, 1000)
}
}

\keyword{ methods }
\keyword{ optimize }
