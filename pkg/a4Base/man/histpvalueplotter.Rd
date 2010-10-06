\name{histpvalueplotter}
\alias{histpvalueplotter}
\title{Workhorse function for the histPvalue function}
\description{
Workhorse function for the histPvalue function. 
This function displays the distribution of the p values using
  a histogram; the horizontal line represents a uniform distribution 
  based on the p value distribution between 0.5 and 1. This represents
  the hypothetical p value distribution arising just by chance.
  This uniform distribution is used to estimate the proportion of differentially
  expressed genes.
}
\usage{
histpvalueplotter(pValue, addLegend = FALSE, xlab = NULL, ylab = NULL, main = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pValue}{numeric vector of p values}
  \item{addLegend}{logical; should a legend be added (TRUE) or not (FALSE; default)}
  \item{xlab}{label for the x axis; defaults to NULL (no label)}
  \item{ylab}{label for the y axis; defaults to NULL (no label)}
  \item{main}{main title for the plot; if NULL (default) no main title is displayed}
  \item{\dots}{further arguments for the \code{hist} call; currently none are used}
}

\author{Willem Talloen and Tobias Verbeke}

\seealso{\code{\link{histPvalue}}, \code{\link{propdegenescalculation}}}
\examples{
if (require(ALL)){
  data(ALL, package = "ALL")
  ALL <- addGeneInfo(ALL)
  ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
 
  tTestResult <- tTest(ALL, "BTtype")
  histPvalue(tTestResult[,"p"], addLegend = TRUE, xlab = "Adjusted P Value")
  histPvalue(tTestResult[,"p"], addLegend = TRUE, main = "Histogram of Adjusted P Values")
  propDEgenesRes <- propDEgenes(tTestResult[,"p"])
}
}
