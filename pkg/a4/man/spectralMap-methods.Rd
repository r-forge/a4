\name{spectralMap-methods}
\docType{methods}
\alias{spectralMap-methods}
\alias{spectralMap,ExpressionSet,character-method}
\title{Methods for Function spectralMap}
\description{
  Methods for spectralMap
}
\section{Methods}{
\describe{

ExpressionSet,character

\item{object = "ExpressionSet", groups = "character", makeLognormal = "logical",
  plot.mpm.args = "list", probe2gene = "logical"}{wrapper around 
  \code{plot.mpm} from the \code{mpm} package}

}}
\arguments{
\item{makeLognormal}{boolean indicating whether one wants to exponentiate the
     data to make them lognormally shaped (\code{TRUE}; the default) or not 
     (\code{FALSE})}
\item{mpm.args}{list of arguments that can be passed to the \code{mpm} function}
\item{plot.mpm.args}{list of arguments that can be passed to the 
     \code{plot.mpm} function that actually draws the plot}
\item{probe2gene}{boolean indicating whether one wants to display the gene symbols
     for the labeled points (\code{TRUE}) or not (\code{FALSE}; the default)}
}
\keyword{methods}
\keyword{hplot}
