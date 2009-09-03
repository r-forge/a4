\name{addQuantilesColors}
\alias{addQuantilesColors}
\title{ Compute quantiles for plotGeneDE function }
\description{
Compute quantiles on mean expression level for plotGeneDE function. Colors of bars in the plot could then be allocated using buckets defined by those quantiles.
}
\usage{
addQuantilesColors(e, ngroups = 3)
}
\arguments{
  \item{e}{ ExpressionSet object to use for computation }
  \item{ngroups}{ Number of groups to be created }
}
\details{
  Number of computed quantiles is equal to (ngroups - 1).
}
\value{
  The ExpressionSet object e is returned, with a new column called colorsQuantilesVector in its slot featureData
}
\author{Eric Lecoutre}
\seealso{ \code{\link{plotLogRatio}}}
\examples{
  data(ALL, package = "ALL")
  ALLQ <- addQuantilesColors(ALL)
  fData(ALLQ)
}
\keyword{manip}

