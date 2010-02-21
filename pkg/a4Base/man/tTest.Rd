\name{tTest}
\alias{tTest}
\title{Use t Test to Compare Two Groups}
\description{
  Use a (modified) t test to compare two groups
}
\usage{
tTest(object, groups, probe2gene = TRUE)
}
\arguments{
  \item{object}{ExpressionSet object}
  \item{groups}{string indicating the name of the variable of the phenoData
    containing the group information}
  \item{probe2gene}{logical; if \code{TRUE} Affymetrix probeset IDs are translated
    into gene symbols; if \code{FALSE} no such translation is conducted}
}
\details{
  For multiple testing the \code{mt.rawp2adjp} function of package
  \code{multtest} is used.
}
\value{
  Object of class \code{"tTest"}, a data frame with the following columns
  \item{gSymbol}{Gene Symbol}
  \item{p}{TODO}
  \item{logRatio}{TODO}
  \item{pBH}{TODO}
  \item{tStat}{TODO}
}
\author{Willem Talloen, Tobias Verbeke}
\examples{
if (require(ALL)){
  data(ALL, package = "ALL")
  ALL <- addGeneInfo(ALL)
  ALL$BTtype <- as.factor(substr(ALL$BT,0,1))

  tTestRes <- tTest(object = ALL,	groups = "BTtype", probe2gene = TRUE)
  volcanoPlot(tTestRes)  
}
}
\seealso{\code{\link[genefilter]{rowttests}}}
\keyword{htest}
