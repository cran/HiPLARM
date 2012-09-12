\name{hiplarSet}
%\docType{methods}
\alias{hiplarSet}
\alias{hiplarSet-methods}
\alias{hiplarSet,character,numeric-method}
\title{ ~~ Methods for Function \code{hiplarSet} in Package \pkg{HiPLARM} ~~}
\description{
 \code{hiplarSet} allows the user to access internal variables within HiPLAR, these
 are generally set automatically but should the user wish they can over ride
 these settings by following the instructions below. The \code{hiplarSet} function
 is useful if users want to choose either the MAGMA or PLASMA libraries or if
 they wish to set the crossover point within the function. The crossover point
 (for users with GPU and CPU support) is the matrix size at which the function switches
 between the PLASMA, CPU libraries and the MAGMA, GPU libraries. This is generally automatic
 but the user may change if they so wish.
}

\section{Methods}{
\describe{


\item{\code{signature(var = "character", val = "numeric")}}{}

}}

\details{
	The \code{var} argument accesses the particular settings that the user wishes to change
	In setting the xover_<function_name> the \code{val} argument should be
	set to a suitable value. These values are generally set automatically during
	setup so the user should be wary when changing it. When accessing
	the library setting the user should choose values from 1 to 3. 1
	sets the PLASMA library for use, 2 sets the MAGMA library and 
	3 sets it to automatic.
}
\examples{
	
		# Sets the PLASMA libary to be used exclusively #
		hiplarSet("hiplar_library",  1)
		# Sets the MAGMA library to be used exclusively #
		hiplarSet("hiplar_library", 2)
		# Enables autotune which selects PLASMA or MAGMA depending on the problem size.#
		hiplarSet("hiplar_library", 3)

		## Methods for setting crossover values ##
		optsize <- 512
		hiplarSet("xover_dgeMatrix_LU", optsize)
		hiplarSet("xover_dgeMatrix_crossprod", optsize)
		hiplarSet("xover_dgeMatrix_dgeMatrix_crossprod", optsize)
		hiplarSet("xover_dgeMatrix_matrix_crossprod", optsize)
		hiplarSet("xover_dgeMatrix_determinant", optsize)
		hiplarSet("xover_dgeMatrix_matrix_mm", optsize)
		hiplarSet("xover_dgeMatrix_matrix_mm",optsize)
		hiplarSet("xover_dgeMatrix_norm",optsize)
		hiplarSet("xover_dgeMatrix_solve", optsize)
		hiplarSet("xover_dgeMatrix_matrix_solve", optsize)
		hiplarSet("xover_dgeMatrix_rcond", optsize)
		hiplarSet("xover_dgeMatrix_LU", optsize)
		hiplarSet("xover_dgeMatrix_crossprod", optsize)
		hiplarSet("xover_dgeMatrix_dgeMatrix_crossprod", optsize)
		hiplarSet("xover_dgeMatrix_matrix_crossprod", optsize)
		hiplarSet("xover_dgeMatrix_determinant", optsize)
		hiplarSet("xover_dgeMatrix_matrix_mm", optsize)
		hiplarSet("xover_dgeMatrix_matrix_mm",optsize)
		hiplarSet("xover_dgeMatrix_norm",optsize)
		hiplarSet("xover_dgeMatrix_solve", optsize)
		hiplarSet("xover_dgeMatrix_matrix_solve", optsize)
		hiplarSet("xover_dgeMatrix_rcond", optsize)
		hiplarSet("xover_dpoMatrix_chol",optsize)
		hiplarSet("xover_dpoMatrix_rcond", optsize)
		hiplarSet("xover_dpoMatrix_solve", optsize)
		hiplarSet("xover_dpoMatrix_dgeMatrix_solve", optsize)
		hiplarSet("xover_dpoMatrix_matrix_solve", optsize)
		hiplarSet("xover_dtrMatrix_chol2inv", optsize)
		hiplarSet("xover_dtrMatrix_dtrMatrix_mm", optsize)
		hiplarSet("xover_dtrMatrix_matrix_mm", optsize)
		hiplarSet("xover_dtrMatrix_solve", optsize)
		hiplarSet("xover_dtrMatrix_matrix_solve", optsize)
		hiplarSet("xover_dsyMatrix_matrix_mm", 0)
		hiplarSet("xover_dsyMatrix_norm", 0)
}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
