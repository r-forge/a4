\name{replicates}
\alias{replicates}
\title{ computes replicates across a vector }
\description{
Given a vector, returns the replicates in order
}
\usage{
replicates(x)
}
\arguments{
  \item{x}{ character or numeric vector }
}
\value{
  numeric vector 
}
\references{R-help mailing list }
\author{Henrique Dallazuanna}
\seealso{ \code{\link{rle}}}
\examples{
  x <- c('a','b','a','a','b','a','c','c','c')
  data.frame(val=x,rep=replicates(x))
}
\keyword{manip}
