\name{tcrossprod}
\title{tcrossproduct using GPU or multi-core CPU}
\alias{tcrossprod}
\alias{tcrossprod,dgeMatrix-method}
\alias{tcrossprod,dtrMatrix-method}
\description{Given matrices \code{x} and \code{y} as arguments, return a matrix
cross-product. 
This is formally equivalent to (but usually slightly faster than) the call \code{t(x) \%*\% y} (tcrossprod) or \code{x \%*\% t(y)} (tcrossprod).
}

\usage{
tcrossprod(x, y)
\S4method{tcrossprod}{dgeMatrix}(x, y)
\S4method{tcrossprod}{dtrMatrix}(x, y)

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
   	\item{tcrossprod}{\code{(signature(x = "dgeMatrix", y = "missing")}: Calls CUBLAS function \code{cublasDsyrk} for GPU enabled
	systems and \code{PLASMA_dsyrk} for multi-core systems. Library settings also affect choice between GPU and CPU. If
	\code{y = NULL}, then it is taken to be the same matrix as \code{x}.}
   
   	\item{tcrossprod}{\code{(signature(x = "dgeMatrix", y = "dgeMatrix")}: Calls MAGMA function \code{magma_dgemm} for GPU enabled
	systems and \code{PLASMA_dgemm} for multi-core systems. Library settings also affect choice between GPU and CPU.}
    
   	\item{tcrossprod}{\code{(signature(x = "dgeMatrix", y = "Matrix")}: Calls MAGMA function \code{magma_dgemm} for GPU enabled
	systems and \code{PLASMA_dgemm} for multi-core systems. Library settings also affect choice between GPU and CPU.}
   
   	\item{tcrossprod}{\code{(signature(x = "dgeMatrix", y = "numeric")}: Calls MAGMA function \code{magma_dgemm} for GPU enabled
	systems and \code{PLASMA_dgemm} for multi-core systems. Library settings also affect choice between GPU and CPU. \code{y} is 
	coerced to a \code{base::matrix}.}
   
   	\item{tcrossprod}{\code{(signature(x = "dgeMatrix", y = "matrix")}: Calls MAGMA function \code{magma_dgemm} for GPU enabled
	systems and \code{PLASMA_dgemm} for multi-core systems. Library settings also affect choice between GPU and CPU.}
	
   	\item{tcrossprod}{\code{(signature(x = "dtrMatrix", y = "dtrMatrix")}: Calls CUBLAS function \code{cublasDtrmm} for GPU enabled
	systems and \code{PLASMA_dtrmm} for multi-core systems. Library settings also affect choice between GPU and CPU.}
   	
	\item{tcrossprod}{\code{(signature(x = "denseMatrix", y = "dtrMatrix")}: \code{y} inherits from virtual class \code{ddenseMatrix}
	Calls MAGMA function \code{magma_dgemm} for GPU enabled systems and \code{PLASMA_dgemm} for multi-core systems. Library settings also affect choice between GPU and CPU.}
   	
\item{tcrossprod}{\code{(signature(x = "matrix", y = "dtrMatrix")}: Calls CUBLAS function \code{cublasDtrmm} for GPU enabled
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
