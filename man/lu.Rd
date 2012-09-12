\name{lu}
\title{(Generalized) Triangular Decomposition of a Matrix}
\alias{lu}
\alias{lu,dgeMatrix-method}
\usage{
lu(x, \dots)
\S4method{lu}{dgeMatrix}(x, warnSing = TRUE, \dots)
}
\description{
  Computes (generalized) triangular decompositions of square and other
  dense matrices using the MAGMA GPU library or the PLASMA library for multi-core CPUs.
  	
}
\arguments{
  \item{x}{a dense matrix.  No missing values or IEEE special values are allowed.}
  \item{warnSing}{(when \code{x} is a
    \code{denseMatrix} logical specifying if a \code{warning}
    should be signalled when \code{x} is singular.}
   \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  \code{lu()} is a generic function with special methods for different types
  of matrices.  Use \code{\link{showMethods}("lu")} to list all the methods
  for the \code{\link{lu}} generic.

  The method for class \code{dgeMatrix} (and all dense
  matrices) is based on the MAGMA \code{"magma_dgetrf"} subroutine and PLASMA
  \code{"PLASMA_dgetrf"}.  It returns
  a decomposition also for singular and non-square matrices. For further 
  details on classes etc. see the Matrix package documentation.
}

\references{
Martin Maechler, Douglas Bates (Matrix package)
}

\examples{

##--- Dense  -------------------------
x <- Matrix(rnorm(9), 3, 3)
lu(x)
dim(x2 <- round(10 * x[,-3]))# non-square
expand(lu2 <- lu(x2))

}
\keyword{array}
\keyword{algebra}
