\name{nlcvTT}
\alias{nlcvTT}
\docType{data}
\title{Data to Demonstrate nlcv and Co Functions}
\description{
  Simulated data set used to demonstrate nlcv and accompanying
  plot functions to study classification problems
}
\usage{data(nlcvTT)}
\format{
  The object is of class \code{"nlcv"}, an object as produced
  by the \code{nlcv} function.
}
\source{
  data simulated using \code{nlcvTT <- nlcv(selBcrAblOrNeg, classVar = 'mol.biol', 
    classdist = "unbalanced", nRuns = 10, fsMethod = "t.test", verbose = TRUE)}
}
\seealso{\code{\link[nlcv]{nlcv}}}
\examples{
  data(nlcvTT)
  if (require(nlcv))
    scoresPlot(nlcvTT, tech = 'svm', nfeat = 25)
}
\keyword{datasets}
