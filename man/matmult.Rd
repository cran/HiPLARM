\name{\%*\%}
\title{Matrix Multiplication of two matrices using GPU or multi-core architectures}
\alias{\%*\%}
\alias{\%*\%,dgeMatrix-method,missing-method}
\alias{\%*\%,dgeMatrix-method,ddenseMatrix-method}
\alias{\%*\%,dgeMatrix-method,matrix-method}
\alias{\%*\%,dpoMatrix-method,missing-method}
\alias{\%*\%,dpoMatrix-method,dgeMatrix-method}
\alias{\%*\%,dpoMatrix-method,matrix-method}
\alias{\%*\%,dtpMatrix-method,ddenseMatrix-method}
\alias{\%*\%,ddenseMatrix-ddenseMatrix}
\alias{\%*\%,ddenseMatrix,dsyMatrix-method}
\alias{\%*\%,ddenseMatrix,dtrMatrix-method}
\alias{\%*\%,dgeMatrix,dgeMatrix-method}
\alias{\%*\%,dgeMatrix,dtpMatrix-method}
\alias{\%*\%,dgeMatrix,matrix-method}
\alias{\%*\%,dMatrix,nMatrix-method}
\alias{\%*\%,dspMatrix,ddenseMatrix-method}
\alias{\%*\%,dspMatrix,matrix-method}
\alias{\%*\%,dsyMatrix,ddenseMatrix-method}
\alias{\%*\%,dsyMatrix,dsyMatrix-method}
\alias{\%*\%,dsyMatrix,matrix-method}
\alias{\%*\%,dtpMatrix,ddenseMatrix-method}
\alias{\%*\%,dtpMatrix,matrix-method}
\alias{\%*\%,dtrMatrix,ddenseMatrix-method}
\alias{\%*\%,dtrMatrix,dtrMatrix-method}
\alias{\%*\%,dtrMatrix,matrix-method}
\alias{\%*\%,Matrix,ANY-method}
\alias{\%*\%,matrix,dgeMatrix-method}
\alias{\%*\%,matrix,dsyMatrix-method}
\alias{\%*\%,matrix,dtpMatrix-method}
\alias{\%*\%,matrix,dtrMatrix-method}
\alias{\%*\%,matrix,Matrix-method}
\alias{\%*\%,Matrix,matrix-method}
\alias{\%*\%,Matrix,numeric-method}
\alias{\%*\%,matrix,pMatrix-method}
\alias{\%*\%,Matrix,pMatrix-method}

\alias{\%*\%,ddenseMatrix-method,ddenseMatrix-method}
\alias{\%*\%,ddenseMatrix-class,ddenseMatrix-class}
\alias{\%*\%,ddenseMatrix,ddenseMatrix-method}
\alias{\%*\%,dtpMatrix-method,matrix-method}
\alias{\%*\%,ddenseMatrix,dsyMatrix-method}
\alias{\%*\%,ddenseMatrix,dtrMatrix-method}
\alias{\%*\%,dtrMatrix-method}
\alias{\%*\%,dtpMatrix-method}	
\alias{\%*\%,dpoMatrix-method}
\alias{\%*\%,dgeMatrix,dgeMatrix-method}
\alias{\%*\%,ddenseMatrix-method}
\alias{\%*\%,dsyMatrix-method}
\alias{\%*\%,dgeMatrix-method}
\description{
	Matrix multiplication using a GPU or multi-core CPU for most dense matrix types
}
\usage{

\S4method{\%*\%}{dpoMatrix}(x, y)
\S4method{\%*\%}{dgeMatrix}(x, y)
\S4method{\%*\%}{dtrMatrix}(x, y)
\S4method{\%*\%}{dtpMatrix}(x, y)
}

\arguments{
  \item{x}{A dense matrix inheriting from \code{Matrix} or of type \code{base::matrix}}
  \item{y}{A dense matrix inheriting from \code{dMatrix} or of type \code{base::matrix}}
}

\section{Methods}{
  \describe{
    
    \item{\%*\%}{\code{signature(x = "ddenseMatrix", y = "ddenseMatrix")}:\code{x} and \code{y}
    inherit from \code{ddenseMatrix}, \code{x} is coerced to \code{dgeMatrix}. For multi-core machines
    the \code{PLASMA_dgemm} is used or for GPU enabled machines \code{magma_dgemm} may be used also.}


    \item{\%*\%}{\code{signature(x = "dgeMatrix", y = "dgeMatrix")}:\code{x} and \code{y} are
    of type \code{dgeMatrix}. For multi-core machines the \code{PLASMA_dgemm} is used or for 
    GPU enabled machines \code{magma_dgemm} may be used also.}
   
   
    \item{\%*\%}{\code{signature(x = "dgeMatrix", y = "matrix")}:For multi-core machines the \code{PLASMA_dgemm} is used or for 
    GPU enabled machines \code{magma_dgemm} may be used also.}
   
    \item{\%*\%}{\code{signature(x = "dsyMatrix", y = "Matrix")}: For multi-core machines the \code{PLASMA_dsymm} is used or for 
    GPU enabled machines \code{cublasDsymm} may be used also.}
   
   
    \item{\%*\%}{\code{signature(x = "dsyMatrix", y = "ddenseMatrix")}:For multi-core machines the \code{PLASMA_dsymm} is used or for 
    GPU enabled machines \code{cublasDsymm} may be used also.}

    \item{\%*\%}{\code{signature(x = "dsyMatrix", y = "dsyMatrix")}:For multi-core machines the \code{PLASMA_dsymm} is used or for 
    GPU enabled machines \code{cublasDsymm} may be used also.}

    \item{\%*\%}{\code{signature(x = "ddenseMatrix", y = "dsyMatrix")}:For multi-core machines the \code{PLASMA_dsymm} is used or for 
    GPU enabled machines \code{cublasDsymm} may be used also.}

    \item{\%*\%}{\code{signature(x = "matrix", y = "dsyMatrix")}:For multi-core machines the \code{PLASMA_dsymm} is used or for 
    GPU enabled machines \code{cublasDsymm} may be used also.}

    \item{\%*\%}{\code{signature(x = "dspMatrix", y = "ddenseMatrix")}:For 
    GPU enabled machines \code{cublasDspmv} may be used. The PLASMA library does not support this type so there is no multi-core support for this function}

    \item{\%*\%}{\code{signature(x = "dspMatrix", y = "Matrix")}:For 
    GPU enabled machines \code{cublasDspmv} may be used also. The PLASMA library does not support this type so there is no multi-core support for this function.}
	
    \item{\%*\%}{\code{signature(x = "dtrMatrix", y = "dtrMatrix")}:For multi-core machines the \code{PLASMA_dtrmm} is used or for 
    GPU enabled machines \code{cublasDtrmm} may be used also.}

    \item{\%*\%}{\code{signature(x = "dtrMatrix", y = "ddenseMatrix")}:For multi-core machines the \code{PLASMA_dtrmm} is used or for 
    GPU enabled machines \code{cublasDtrmm} may be used also.}
    
    \item{\%*\%}{\code{signature(x = "dtrMatrix", y = "Matrix")}:For multi-core machines the \code{PLASMA_dtrmm} is used or for 
    GPU enabled machines \code{cublasDtrmm} may be used also.}

    \item{\%*\%}{\code{signature(x = "ddenseMatrix", y = "dtrMatrix")}:For multi-core machines the \code{PLASMA_dtrmm} is used or for 
    GPU enabled machines \code{cublasDtrmm} may be used also.}
    
    \item{\%*\%}{\code{signature(x = "matrix", y = "dtrMatrix")}:For multi-core machines the \code{PLASMA_dtrmm} is used or for 
    GPU enabled machines \code{cublasDtrmm} may be used also.}
	
    \item{\%*\%}{\code{signature(x = "dtpMatrix", y = "ddenseMatrix")}:For 
    GPU enabled machines \code{cublasDtpmv} may be used. The PLASMA library does not support this type so there is no multi-core support for this function}

    \item{\%*\%}{\code{signature(x = "dgeMatrix", y = "dtpMatrix")}:For 
    GPU enabled machines \code{cublasDtpmv} may be used also. The PLASMA library does not support this type so there is no multi-core support for this function.}

   }


}

\details{
	For further details on classes and methods see the full
	Matrix package documentation.
}

\references{
Martin Maechler, Douglas Bates (Matrix package)
}

\examples{
 p <- 256
 X <- Matrix(rnorm(p*p), p, p) # random square matrix for large p
 Y <- Matrix(rnorm(p*p), p ,p)
 Z <- X %*% Y

 #dtr triangular Matrix
 X <- triu(X)
 Y <- triu(y)
 Z <- X %*% Y
}
\keyword{algebra}
\keyword{array}
