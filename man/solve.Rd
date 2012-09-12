\name{solve}
\title{Solve a linear system Ax=b using GPU or multi-core architectures}
\alias{solve}
\alias{solve,dgeMatrix-method,missing-method}
\alias{solve,dgeMatrix-method,ddenseMatrix-method}
\alias{solve,dgeMatrix-method,matrix-method}
\alias{solve,dgeMatrix-method,sparseMatrix-method}
\alias{solve,dpoMatrix-method,missing-method}
\alias{solve,dpoMatrix-method,dgeMatrix-method}
\alias{solve,dpoMatrix-method,matrix-method}
\alias{solve,dtpMatrix-method,ddenseMatrix-method}
\alias{solve,dtpMatrix-method,matrix-method}

\alias{solve,dtrMatrix-method}
\alias{solve,dtpMatrix-method}	
\alias{solve,dpoMatrix-method}
\alias{solve,dgeMatrix-method}


\description{
  Solve, using a GPU or multi-core CPU, a linear system Ax=b where A is one of \code{dgeMatrix},
  \code{dpoMatrix}, \code{dtrMatrix} or \code{dtpMatrix}.
}
\usage{

\S4method{solve}{dpoMatrix}(A, b, \dots)
\S4method{solve}{dgeMatrix}(A, b, \dots)
\S4method{solve}{dtrMatrix}(A, b, \dots)
\S4method{solve}{dtpMatrix}(A, b, \dots)

}

\arguments{
  \item{A}{Square dense matrix inheriting from dMatrix.}
  \item{b}{Matrix/vector inheriting from dMatrix or base::matrix.}
  \item{\dots}{potentially further arguments passed to methods.}
}

\section{Methods}{
  \describe{
    
    \item{solve}{\code{signature(A = "dgeMatrix", b = "missing")}: Sets \code{b} as
    the identity matrix and calculates inverse of \code{A}; Uses PLASMA_dgetri or magma_dgetri
    for multi-core CPU and GPU respectively. Also uses some LAPACK \code{dgetri} to test singularity.}

	\item{solve}{\code{signature(A = "dgeMatrix", b = "ddenseMatrix")}: Solves a linear system 
	where b inherits from \code{ddenseMatrix}. If the routine is calling MAGMA \code{"magma_dgetrs_gpu"} is 
	called. If PLASMA is chosen \code{"PLASMA_dgetrs"} is called.}

	\item{solve}{\code{signature(A = "dgeMatrix", b = "matrix")}:  Solves a linear system 
	where b is of type matrix from the R \code{base}. If the routine is calling MAGMA \code{"magma_dgetrs_gpu"} is 
	called. If PLASMA is chosen \code{"PLASMA_dgetrs"} is called.}

	
	\item{solve}{\code{signature(A = "dgeMatrix", b = "sparseMatrix")}:  Solves a linear system 
	where b is a sparseMatrix from the Matrix package. For this routine the sparse matrix is coerced to a real dense matrix. 
	If the routine is calling MAGMA \code{"magma_dgetrs_gpu"} is 
	called. If PLASMA is chosen \code{"PLASMA_dgetrs"} is called.}
   
   	\item{solve}{\code{signature(A = "dpoMatrix", b = "missing")}: Sets \code{b} as
    the identity matrix and calculates inverse of \code{A}; Uses \code{"PLASMA_dgetri"} or \code{"magma_dpotri_gpu"}
    for multi-core CPU and GPU respectively.}

	\item{solve}{\code{signature(A = "dpoMatrix", b = "dgeMatrix")}: Solves a linear system 
	where b is of type \code{dgeMatrix}. If the routine is calling MAGMA \code{"magma_dpotrs_gpu"} is 
	called. If PLASMA is chosen \code{"PLASMA_dpotrs"} is called.}

	\item{solve}{\code{signature(A = "dpoMatrix", b = "matrix")}:  Solves a linear system 
	where b is of type matrix from the R \code{base}. If the routine is calling MAGMA \code{"magma_dpotrs_gpu"} is 
	called. If PLASMA is chosen \code{"PLASMA_dpotrs"} is called.}

 	\item{solve}{\code{signature(A = "dtpMatrix", b = "missing")}: Sets \code{b} as
    the identity matrix and calculates inverse of \code{A}; This is only supported on the GPU
    so there is no PLASMA call here. Also the CUBLAS library is called here using \code{cublasDtpsv}}.

	\item{solve}{\code{signature(A = "dtpMatrix", b = "ddenseMatrix")}: Solves a linear system 
	where b inherits from \code{ddenseMatrix}. Again there is no MAGMA or PLASMA support for the \code{dtpMatrix} type
	but for GPU capable systems we call \code{cublasDtpsv}.}
	
	\item{solve}{\code{signature(A = "dtrMatrix", b = "missing")}: Sets \code{b} as
    the identity matrix and calculates inverse of \code{A}; Uses \code{"PLASMA_dtrtri"} or \code{"magma_dtrtri"}
    for multi-core CPU and GPU respectively.}

	\item{solve}{\code{signature(A = "dtrMatrix", b = "ddenseMatrix")}: Solves a linear system where b inherits from \code{ddenseMatrix}. If the routine is calling MAGMA \code{"magma_dtrsm"} is called. If PLASMA is chosen \code{"PLASMA_dtrsm"} is called.}

	\item{solve}{\code{signature(A = "dtrMatrix", b = "matrix")}:  Solves a linear system 
	where b is of type matrix from the R \code{base}. If the routine is calling MAGMA \code{"magma_dtrsm"} is 
	called. If PLASMA is chosen \code{"PLASMA_dtrsm"} is called.}

	\item{solve}{\code{signature(A = "dtrMatrix", b = "Matrix")}:  Solves a linear system 
	where b inherits from \code{Matrix}. If the routine is calling MAGMA \code{"magma_dtrsm"} is 
	called. If PLASMA is chosen \code{"PLASMA_dtrsm"} is called.}
   
	\item{solve}{\code{signature(A = "dtrMatrix", b = "dMatrix")}:  Solves a linear system 
	where b inherits from \code{Matrix}. If the routine is calling MAGMA \code{"magma_dtrsm"} is 
	called. If PLASMA is chosen \code{"PLASMA_dtrsm"} is called.}

   }
}
\references{
Martin Maechler, Douglas Bates (Matrix package)
}

\details{
	For further details on classes and methods see the full
	Matrix package documentation.
}

\examples{
 p <- 128
 A <- Matrix(rnorm(p*p), p, p) # random square matrix for large p
 x_init <- vector("numeric", p)
 b <- A %*% x_init
 x <- solve(A, b)
 stopifnot(identical(x, x_init))

}
\keyword{algebra}
\keyword{array}
