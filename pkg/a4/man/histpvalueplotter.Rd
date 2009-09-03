\name{histpvalueplotter}
\alias{histpvalueplotter}
%- Also NEED an '\alias' for EACH other topic documented here.
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
histpvalueplotter(pValue, addLegend = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pValue}{ ~~Describe \code{pValue} here~~ }
  \item{addLegend}{ ~~Describe \code{addLegend} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}

\author{Willem Talloen and Tobias Verbeke}

\seealso{\code{\link{histPvalue}}, \code{\link{propdegenescalculation}}}
\examples{
data(ALL, package = "ALL")
ALL <- addGeneInfo(ALL)
ALL$BTtype <- as.factor(substr(ALL$BT,0,1))

tTestResult <- tTest(ALL, "BTtype")
histPvalue(tTestResult[,"p"], addLegend = TRUE)
propDEgenesRes <- propDEgenes(tTestResult[,"p"])
}
%% Add one or more standard keywords, see file 'KEYWORDS' in the
%% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
