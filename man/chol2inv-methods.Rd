\name{chol2inv-methods}
\docType{methods}
\alias{chol2inv}
\alias{chol2inv-methods}
\alias{chol2inv,dtrMatrix-method}

\title{Inverse from Cholesky}
\description{
  Invert a symmetric, positive definite square matrix from its Cholesky
  decomposition.  Equivalently, compute \eqn{(X'X)^{-1}}{(X'X)^(-1)}
  from the (\eqn{R} part) of the QR decomposition of \eqn{X}.
  \cr
  Even more generally, given an upper triangular matrix \eqn{R},
  compute \eqn{(R'R)^{-1}}{(R'R)^(-1)}.
}
\usage{
\S4method{chol2inv}{dtrMatrix}(x, \dots)
}
\arguments{
   \item{x}{a matrix(-like) object; see below.}
   \item{\dots}{not used here; for compatibility with other methods.}
 }
\section{Methods}{
  \describe{

    \item{x = "dtrMatrix"}{method for the numeric triangular matrices,
      built on the MAGMA \command{magma_dpotri} and PLASMA \command{PLASMA_dpotri}.}


  }
}
\references{
Martin Maechler, Douglas Bates (Matrix package)
}

\examples{
(M  <- Matrix(cbind(1, 1:3, c(1,3,7))))
(cM <- chol(M)) # a "Cholesky" object, inheriting from "dtrMatrix"
chol2inv(cM) \%*\% M # the identity
stopifnot(all(chol2inv(cM) \%*\% M - Diagonal(nrow(M))) < 1e-10)
}
\keyword{methods}
\keyword{algebra}

