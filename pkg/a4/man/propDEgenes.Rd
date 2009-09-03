\name{propDEgenes}
\alias{propDEgenes}
\title{Generic function to compute the proportion of differentially expressed
  genes that are present}
\description{
  Generic function to compute the proportion of differentially expressed
  genes that are present. Methods are available for objects of class
  'tTest', 
}
\usage{
propDEgenes(object, ...)
}
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
