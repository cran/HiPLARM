\name{checkFile}
\alias{checkFile}
\title{Startup function that reads in the optimised crossover points}
\description{
	This function is run automatically at load time and is the function
	that reads in the already generated optimised crossover points or a default value.}
\usage{
checkFile()
}

\details{
	This function is run automatically at load time and is the function
	that reads in the already generated optimised crossover points. If these
	do not exist a default value is read in. The user is also informed of how many
	functions have been optimised. 

}
\examples{
	checkFile()
}
