\name{probe2gene}
\alias{probe2gene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Translate Affymetrix probeset IDs into gene symbols}
\description{Translate Affymetrix probeset IDs into gene symbols}
\usage{
probe2gene(probesetIds, chipPkg)
}
\arguments{
  \item{probesetIds}{Affymetrix probeset IDs}
  \item{chipPkg}{string indicating the annotation package for the chip}
}
\value{
  Vector containing the respective gene symbols
}
\author{Tobias Verbeke}

\seealso{\code{\link{spectralMap}}, \code{\link[a4Classif]{lassoClass}}, ...}
\examples{
if (require(ALL)){
  data(ALL, package = "ALL")
  chip <- annotation(ALL)
  chipAnnotationPkg <- paste(chip, "db", sep = ".")
  res <- probe2gene(featureNames(ALL), chipAnnotationPkg)
  head(res)
}
}
\keyword{manip}