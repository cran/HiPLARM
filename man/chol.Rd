\name{chol}
\title{Cholesky Decomposition using GPU or multi-core CPU}
\alias{chol}
%\alias{chol-methods}
%\alias{chol,Matrix-method}
%\alias{chol,ddenseMatrix-method}
%\alias{chol,dgeMatrix-method}
\alias{chol,dpoMatrix-method}
%\alias{chol,dppMatrix-method}

\description{
  Compute the Cholesky factorization of a real symmetric
  positive-definite square matrix using GPU or multi-core CPU.
}
\usage{
chol(x, \dots)
\S4method{chol}{dpoMatrix}(x, pivot = FALSE, \dots)
}
\arguments{
  \item{x}{if \code{x} is not positive definite,
    an error is signalled.}
  \item{pivot}{logical indicating if pivoting is used. This is not supported
  for the GPU implementation}
  \item{\dots}{potentially further arguments passed to methods.}
}
\details{
	For further details on classes and methods see the full
	Matrix package documentation.
}
\section{Methods}{
  \describe{
    \item{chol}{\code{signature(x = "dgeMatrix")}: works via
      \code{"dpoMatrix"}}

    \item{chol}{\code{signature(x = "dpoMatrix")}:
      Returns (and stores) the Cholesky decomposition of \code{x}, via
      LAPACK routines \code{dlacpy} and either \code{magma_dpotrf} or \code{PLASMA_dpotrf}.}

   }
}

\references{
Martin Maechler, Douglas Bates (Matrix package)
}

\examples{
showMethods(chol, inherited = FALSE) # show different methods

sy2 <- new("dsyMatrix", Dim = as.integer(c(2,2)), x = c(14, NA,32,77))
(c2 <- chol(sy2))#-> "Cholesky" matrix
stopifnot(all.equal(c2, chol(as(sy2, "dpoMatrix")), tol= 1e-13))
str(c2)

## An example where chol() can't work
(sy3 <- new("dsyMatrix", Dim = as.integer(c(2,2)), x = c(14, -1, 2, -7)))
try(chol(sy3)) # error, since it is not positive definite
}
\keyword{algebra}
\keyword{array}
