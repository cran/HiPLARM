\name{rcond}
\title{Estimate the Reciprocal Condition Number using GPU and multi-core CPU}
\alias{rcond}
\alias{rcond,dpoMatrix-method}
\alias{rcond,dgeMatrix-method}

\description{
  Estimate the reciprocal of the condition number of a matrix using
  the MAGMA library for GPU or PLASMA library for multi-core CPUs.
}

\usage{
\S4method{rcond}{dgeMatrix}(x, norm, \dots)
\S4method{rcond}{dpoMatrix}(x, norm, \dots)
}

%\section{Methods}{
%\describe{

%\item{\code{signature(x = "dgeMatrix", y = "character")}}{}
%\item{\code{signature(x = "dpoMatrix", y = "character")}}{}
%\item{\code{signature(x = "dpoMatrix", y = "missing")}}{}

%}
%}

\arguments{
  \item{x}{an \R object that inherits from the \code{Matrix} class.}
  \item{norm}{
    Character indicating the type of norm to be used in the estimate.
    The default is \code{"O"} for the 1-norm (\code{"O"} is equivalent
    to \code{"1"}).  The other possible value is \code{"I"} for the
    infinity norm, see also \code{norm}.
  }
  
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{ Uses the LU decomposition using \code{magma_dgetrf} or
\code{PLASMA_dgetrf} with the respective MAGMA or PLASMA norms if
available (see norm documentation for more info).It also calls LAPACK \code{dgecon} 
as this is not supported in MAGMA/PLASMA. For further details on classes and methods 
see the Matrix package documentation.
}
\value{
  An estimate of the reciprocal condition number of \code{x}.
}

\references{
Martin Maechler, Douglas Bates (Matrix package)
}

\examples{
x <- Matrix(rnorm(9), 3, 3)
rcond(x)
## typically "the same" (with more computational effort):
1 / (norm(x) * norm(solve(x)))
rcond(Hilbert(9))  # should be about 9.1e-13

## For non-square matrices:
rcond(x1 <- cbind(1,1:10))# 0.05278
rcond(x2 <- cbind(x1, 2:11))# practically 0, since x2 does not have full rank

}
\keyword{array}
\keyword{algebra}
