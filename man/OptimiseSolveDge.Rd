\name{OptimiseSolveDge}
\alias{OptimiseSolveDge}
\title{Optimise the solve routine for a dge matrix}
\description{
	This computes the optimal crossover point 
	between the PLASMA (CPU) and MAGMA (GPU)
	libraries. For smaller matrices PLASMA is generally used
	and for larger matrices MAGMA will be used. This function
	calculates this for solve function in this package and saves
	it for future use.
}
\usage{
OptimiseSolveDge(increment = 128, verbose = FALSE)
}

\arguments{
\item{increment}{This is the step between the problem sizes being calculate. We compare the timings at different
problem sizes for the PLASMA and MAGMA implementations. This is the size between those problem sizes. A smaller increment
will give a more accurate crossover point, however, 128 should be sufficient.}
\item{verbose}{This displays the timing information for the operations if it is set to TRUE. It is FALSE by default.}
}
\details{
	This simply runs the solve function, for a dge matrix, testing it numerous times at different problem sizes. It benchmarks the PLASMA 
	and MAGMA libraries against each other. For smaller problem sizes PLASMA will be faster, as the size increases the
	MAGMA library will become more optimal. This routine will find and save that optimal point for future use.
}
\examples{
	OptimiseSolveDge(256,TRUE)
}
