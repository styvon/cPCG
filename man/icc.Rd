\name{icc}
\alias{icc}

\title{
Incomplete Cholesky Factorization
}
\description{
Incomplete Cholesky factorization method to generate preconditioning matrix for conjugate gradient method.  
}
\usage{
icc(A)
}

\arguments{
  \item{A}{
  matrix, symmetric and positive definite.
  }
}
\details{
Performs incomplete Cholesky factorization on the input matrix A, the output matrix is used for preconditioning in pcgsolve() if "ICC" is specified as the preconditioner.    
}


\value{
Returns a matrix after incomplete Cholesky factorization.
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
out <- icc(test_A)
}
}

\keyword{ methods }
