\name{determinant}
\title{Calculate the determinant using GPU and multi-core CPU}
\alias{determinant,dgeMatrix-method,missing-method}
\alias{determinant,dgeMatrix-method,logical-method}
\alias{determinant,dgeMatrix-method}
\alias{determinant,dgeMatrix-method}

\description{
  Estimate the determinant of a matrix using
  the MAGMA library for GPU or PLASMA library for multi-core CPUs.
}

\usage{
\S4method{determinant}{dgeMatrix}(x, logarithm, \dots)
\S4method{determinant}{dgeMatrix}(x, logarithm, \dots)
}


\section{Methods}{
\describe{
	\item{determinant}{\code{signature(x = "Matrix", logarithm = "missing")}: and}
	\item{determinant}{\code{signature(x = "Matrix", logarithm = "logical")}:compute the (\eqn{\log}) determinant of \code{x}.  
	The method chosen depends on the actual Matrix class of \code{x}.  
	Note that \code{base::det} also works for all our matrices, calling the appropriate \code{determinant()} method.  
	The \code{Matrix::det} is an exact copy of \code{base::det}, but in the correct namespace, and hence calling the S4-aware version of
	\code{determinant()}.}
}
}

\arguments{
  \item{x}{an \R object of type \code{dgeMatrix}.}
  \item{logarithm}{logical; if \code{TRUE} (default) return the logarithm of the modulus of the determinant.
  }
  
  \item{\dots}{No further arguments at present.}
}
\details{ Uses the LU decomposition using \code{magma_dgetrf} or
\code{PLASMA_dgetrf}.For further details
and methods see the Matrix package documentation
or indeed the base package.
}

\references{
Martin Maechler, Douglas Bates (Matrix package)
}

\examples{
	slotNames("Matrix")

	cl <- getClass("Matrix")
	names(cl@subclasses) # more than 40 ..

	showClass("Matrix")#> output with slots and all subclasses

	(M <- Matrix(c(0,1,0,0), 6, 4))
	dim(M)
	diag(M)
	cm <- M[1:4,] + 10*Diagonal(4)
	diff(M)
## can reshape it even :
	dim(M) <- c(2, 12)
	M
	stopifnot(identical(M, Matrix(c(0,1,0,0), 2,12)),all.equal(det(cm),
	determinant(as(cm,"matrix"), log=FALSE)$modulus, check.attr=FALSE))
}
\keyword{array}
\keyword{algebra}
