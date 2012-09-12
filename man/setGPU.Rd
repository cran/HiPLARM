\name{setGPU}
\docType{methods}
\alias{setGPU}
\alias{setGPU,logical,character-method}
\alias{setGPU,logical,missing-method}
\alias{setGPU,numeric,character-method}
\alias{setGPU,numeric,missing-method}
\title{ ~~ Methods for Function \code{setGPU} in Package \pkg{HiPLARM} ~~}
\description{
 ~~ Methods for function \code{setGPU} in package \pkg{HiPLARM} ~~
}
\usage{
setGPU(proc, interface)

}
\arguments{
	\item{proc}{A numeric or logical value that enables
	or disable the argument passed to \code{interface}, i.e the CPU or GPU interface
	or if no argument is passed to interface then the MAGMA library is enabled or disabled.}
	\item{interface}{A character string of either "CPU" or "GPU" determining which
	interface to use. If no argument is given the MAGMA library is disabled.}
}


\details{
	The MAGMA library provides two different interfaces to the GPU. The 
	"CPU" interface automatically transfers the data to and from the GPU
	as well as choosing the multi-gpu option automatically should the user have multiple
	GPUs. This is the most flexible or options and should it be available is
	used by default. However, some MAGMA functions are only implemented for
	the "GPU" interface. In this case the data transfers and other options must be
	given explicitly for the function call. If this is the case, this is the default option, 
	however, at an R level the user will see no difference. 
	We also provide the option of disabling the GPU and 
	using the regular LAPACK call should the user wish to test results.
}
\section{Methods}{
\describe{

\item{\code{signature(proc = "logical", interface = "character")}}{
	This allows the user to choose the CPU or GPU interface. Should
	the user wish to enable a certain interface the string
	"CPU" or "GPU" are passed as the \code{interface} argument and
	the \code{proc} argument is set to \code{TRUE} or \code{FALSE}.
}

\item{\code{signature(proc = "logical", interface = "missing")}}{
	This allows the user to enable or disable the GPU if
	no argument is passed to \code{interface}. Setting \code{proc}
	to \code{TRUE} or \code{FALSE} enables or disables the GPU.
}

\item{\code{signature(proc = "numeric", interface = "character")}}{
	This allows the user to choose the CPU or GPU interface. Should
	the user wish to enable a certain interface the string
	"CPU" or "GPU" are passed as the \code{interface} argument and
	the \code{proc} argument is set to \code{1} or \code{0}.
}

\item{\code{signature(proc = "numeric", interface = "missing")}}{
	This allows the user to enable or disable the GPU if
	no argument is passed to \code{interface}. Setting \code{proc}
	to \code{1} or \code{0} enables or disables the GPU.
}
}}
\keyword{methods}

\examples{
setGPU(TRUE,"CPU") # enables CPU interface where available
setGPU(TRUE,"GPU") # enables GPU interface only
setGPU(FALSE) #disables MAGMA library
}
