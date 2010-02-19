\name{propDEgenes-methods}
\alias{propDEgenes-methods}
\alias{propDEgenes,limma-method}
\alias{propDEgenes,numeric-method}
\title{Methods for propDEgenes}
\description{
  Methods for propDEgenes 
}
\section{Methods}{
\describe{

limma

\item{object = "limma"}{propDEgenes method for a limma object}

numeric

\item{object = "numeric"}{propDEgenes method for a numeric vector, i.e. a vector
  of P Values}

}}
\arguments{
  \item{object}{object of class}
  \item{\dots}{further arguments for the method (currently none implemented)}
}
\value{
  numeric of length one giving the proportion of differentially expressed
  genes
}
% TODO add reference book ?
\author{Willem Talloen and Tobias Verbeke}
\seealso{\code{\link{propDEgenes-methods}}}
% TODO add examples section
\keyword{htest}
