\name{norm}
\title{Matrix Norms}
\alias{norm}
\alias{norm,dMatrix-method}
\alias{norm,dgeMatrix-method}
\alias{norm,dpoMatrix-method}
\alias{norm,dsyMatrix-method}

\alias{norm,ldenseMatrix-method}
\alias{norm,ndenseMatrix-method}
\alias{norm,ANY,missing-method}
\alias{norm,matrix,character-method}
\alias{norm,Matrix,character-method}

\description{
  Computes a matrix norm of \code{x}, using MAGMA for GPUs or PLASMA
  for multi-core CPUs.
  The norm can be the one norm, the infinity norm, the Frobenius norm,
  or the maximum modulus among elements of a matrix, as determined by
  the value of \code{type}. Not all norms are supported by MAGMA or PLASMA and 
  these are documented below. 
}
\usage{
\S4method{norm}{dsyMatrix}(x, type, \dots)
\S4method{norm}{dgeMatrix}(x, type, \dots)
\S4method{norm}{dpoMatrix}(x, type, \dots)
\S4method{norm}{ldenseMatrix}(x, type, \dots)
\S4method{norm}{ndenseMatrix}(x, type, \dots)
}
\arguments{
  \item{x}{
    a real matrix.
  }
  \item{type}{
    A character indicating the type of norm desired. All these norms 
    are supported in the latest PLASMA library (v2.4.6 asof writing).

    \describe{
      \item{\code{"O"}, \code{"o"} or \code{"1"}}{specifies the one norm,
	(maximum absolute column sum); This is not supported in MAGMA currently.}
      \item{\code{"I"} or \code{"i"}}{specifies the infinity norm (maximum
	absolute row sum);}
      \item{\code{"F"} or \code{"f"}}{specifies the Frobenius norm (the
	Euclidean norm of \code{x} treated as if it were a vector); This is not
	supported in PLASMA versions < 2.4.6 or in MAGMA.}
      \item{\code{"M"} or \code{"m"}}{specifies the maximum modulus of
	all the elements in \code{x}. This is not supported in MAGMA for
	dgeMatrix, but is for symmetric matrices \code{dpoMatrix} and \code{dsyMatrix}.}
    }
    The default is \code{"O"}.  Only the first character of
    \code{type[1]} is used.
  }
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
    A numeric value of class \code{"norm"}, representing the quantity
    chosen according to \code{type}.
}
\details{
  When or if the GPU is used \code{magmablas_dlange} for \code{dgeMatrix} is called
  or \code{magmablas_dlansy} is used for symmetric \code{dsyMatrix} or \code{dpoMatrix}.
  When the multi-core library PLASMA is used \code{PLASMA_dlange} and \code{PLASMA_dlansy}
  are used for the respective matrix types mentioned previously.
}

\references{
Martin Maechler, Douglas Bates (Matrix package) 
}

\examples{
x <- Hilbert(9)
norm(x, "1")
norm(x, "I")
norm(x, "F")
norm(x, "M")
}
\keyword{algebra}
