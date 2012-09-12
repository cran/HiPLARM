\name{crossprod}
\title{Crossproduct using GPU or multi-core CPU}
\alias{crossprod}
\alias{crossprod,dgeMatrix-method}
\alias{crossprod,dtrMatrix-method}
\description{Given matrices \code{x} and \code{y} as arguments, return a matrix
cross-product. 
This is formally equivalent to (but usually slightly faster than) the call \code{t(x) \%*\% y} (crossprod) or \code{x \%*\% t(y)} (tcrossprod).
}

\usage{
crossprod(x, y)
\S4method{crossprod}{dgeMatrix}(x, y)
\S4method{crossprod}{dtrMatrix}(x, y)

}
\arguments{
  \item{x}{dense matrix or vector represented as a single row matrix.}
  \item{y}{dense matrix or vector represented as single row matrix.}
}
\details{
	For further details on classes and methods see the full
	Matrix package documentation.
}
\section{Methods}{
  \describe{
   	\item{crossprod}{\code{(signature(x = "dgeMatrix", y = "missing")}: Calls CUBLAS function \code{cublasDsyrk} for GPU enabled
	systems and \code{PLASMA_dsyrk} for multi-core systems. Library settings also affect choice between GPU and CPU. If
	\code{y = NULL}, then it is taken to be the same matrix as \code{x}.}
   
   	\item{crossprod}{\code{(signature(x = "dgeMatrix", y = "dgeMatrix")}: Calls MAGMA function \code{magma_dgemm} for GPU enabled
	systems and \code{PLASMA_dgemm} for multi-core systems. Library settings also affect choice between GPU and CPU.}
    
   	\item{crossprod}{\code{(signature(x = "dgeMatrix", y = "Matrix")}: Calls MAGMA function \code{magma_dgemm} for GPU enabled
	systems and \code{PLASMA_dgemm} for multi-core systems. Library settings also affect choice between GPU and CPU.}
   
   	\item{crossprod}{\code{(signature(x = "dgeMatrix", y = "numeric")}: Calls MAGMA function \code{magma_dgemm} for GPU enabled
	systems and \code{PLASMA_dgemm} for multi-core systems. Library settings also affect choice between GPU and CPU. \code{y} is 
	coerced to a \code{base::matrix}}
   
   	\item{crossprod}{\code{(signature(x = "dgeMatrix", y = "matrix")}: Calls MAGMA function \code{magma_dgemm} for GPU enabled
	systems and \code{PLASMA_dgemm} for multi-core systems. Library settings also affect choice between GPU and CPU.}
	
   	\item{crossprod}{\code{(signature(x = "dtrMatrix", y = "dtrMatrix")}: Calls the CUBLAS function \code{cublasDtrmm} for GPU enabled
	systems and \code{PLASMA_dtrmm} for multi-core systems. Library settings also affect choice between GPU and CPU.}
   	
	\item{crossprod}{\code{(signature(x = "dtrMatrix", y = "ddenseMatrix")}: \code{y} inherits from virtual class \code{ddenseMatrix}
	Calls CUBLAS function \code{magma_dgemm} for GPU enabled systems and \code{PLASMA_dgemm} for multi-core systems. Library settings also affect choice between GPU and CPU.}
   	
\item{crossprod}{\code{(signature(x = "dtrMatrix", y = "matrix")}: Calls CUBLAS function \code{cublasDtrmm} for GPU enabled
systems and \code{PLASMA_dtrmm} for multi-core systems. Library settings also affect choice between GPU and CPU.}
   }
}
\references{
Martin Maechler, Douglas Bates (Matrix package)
}
\examples{
m <- matrix(0, 400, 500)
set.seed(12)
m[runif(314, 0, length(m))] <- 1
mm <- as(m, "dgeMatrix")
object.size(m) / object.size(mm) # smaller by a factor of > 200

## tcrossprod() is very fast:
system.time(tCmm <- tcrossprod(mm))
system.time(cm <- crossprod(t(m)))
system.time(cm. <- tcrossprod(m))

stopifnot(cm == as(tCmm, "matrix"))

}
\keyword{algebra}
\keyword{array}
