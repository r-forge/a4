\name{histPvalue}
\alias{histPvalue}
\alias{histPvalue,tTest-method}
\alias{histPvalue,limma-method}
\alias{histPvalue,MArrayLM-method}
\alias{histPvalue,numeric-method}
\title{Plot the Distribution of P Values}
\description{
  This function displays the distribution of the p values using
  a histogram; the horizontal line represents a uniform distribution 
  based on the p value distribution between 0.5 and 1. This represents
  the hypothetical p value distribution arising just by chance.
  This uniform distribution is used to estimate the proportion of differentially
  expressed genes.
}
\usage{
histPvalue(object,  ...)
\S4method{histPvalue}{MArrayLM}(object, coef, ...)
}
\arguments{
  \item{object}{either a numeric vector of p-values,
  or an object of class \code{tTest}, \code{limma} or \code{MArrayLM}}
  \item{coef}{index of the coefficient for which the p values should be plotted;
    only applies to the MArrayLM method}
  \item{\dots}{further arguments passed to the method}
}
\value{
  The histogram is displayed on the current device. 
}
\references{
  Goehlmann, H. and W. Talloen (2009). Gene Expression Studies Using Affymetrix
     Microarrays, Chapman \& Hall/CRC, p. 253.
}
\author{Willem Talloen and Tobias Verbeke}
\seealso{\code{\link{hist}}, \code{\link{histpvalueplotter}}}
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
\keyword{dplot}
