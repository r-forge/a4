\name{propdegenescalculation}
\alias{propdegenescalculation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of proportion of differentially expressed genes}
\description{
Estimation of proportion of differentially expressed genes.
 This estimation is based on a histogram of the p-values. More specifically,
  based on the horizontal line representing a uniform distribution 
  based on the p value distribution between 0.5 and 1. This represents
  the hypothetical p value distribution arising just by chance.
 All genes with small p-values above this line reflect
  the expected number of differentially expressed genes not by chance. }
\usage{
propdegenescalculation(pValue)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pValue}{a vector of p-values}
}

\author{Willem Talloen and Tobias Verbeke}

\seealso{\code{\link{histPvalue}}}
\examples{
if (require(ALL)){
  data(ALL, package = "ALL")
  ALL <- addGeneInfo(ALL)
  ALL$BTtype <- as.factor(substr(ALL$BT,0,1))

  tTestResult <- tTest(ALL, "BTtype")
  histPvalue(tTestResult[,"p"], addLegend = TRUE)
  propDEgenesRes <- propDEgenes(tTestResult[,"p"])
}
}

