\name{HiPLARM}
\alias{HiPLARM}
\alias{HiPLARM}
\docType{package}
\title{
	HiPLARM
}
\description{
	Multi-core and GPU support for linear algebra functions in the recommended R package Matrix.
}
\details{
\tabular{ll}{
Package: \tab HiPLARM\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2012-08-04\cr
License: \tab GPL(>=2)\cr
}

The HiPLARM package can be used in exactly the same manner as the Matrix package on which it depends.
The functionality maps that of the Matrix package exactly bar some specialised functions which are documented here. These Optimise functions
are only run once to allow HiPLARM to be optimised for particular systems.
}

\author{
Peter Nash, Vendel Szeremi

Maintainer: Giovanni Montana  <support@hiplar.org>

}

\references{
\url{http://icl.cs.utk.edu/plasma}

\url{http://icl.cs.utk.edu/magma/}

Song, F., Tomov, S., Dongarra, J. (2012)
\emph{Enabling and Scaling Matrix Computations on Heterogeneous Multi-Core and Multi-GPU Systems,}
\dQuote{26th ACM International Conference on Supercomputing (ICS 2012), ACM, San Servolo Island, Venice, Italy}

Kurzak, J., Luszczek, P., Faverge, M., Dongarra, J. (April 2012)
\emph{Programming the LU Factorization for a Multicore System with Accelerators}
\dQuote{Proceedings of VECPAR 12, Kobe, Japan}

Haidar, A., Ltaief, H., Dongarra, J. (2011)
\emph{Parallel Reduction to Condensed Forms for Symmetric Eigenvalue Problems using Aggregated Fine-Grained and Memory-Aware Kernels,}
\dQuote{Proceedings of 2011 International Conference for High Performance Computing, Networking, Storage and Analysis (SC11), Seattle, WA}
}
\keyword{ package }

\examples{
library(HiPLARM)
x <- Matrix(rnorm(3 * 3), ncol = 3)
y <- Matrix(rnorm(3 * 3), ncol = 3)
z <- x\%*\%y
y <- lu(z)
}
