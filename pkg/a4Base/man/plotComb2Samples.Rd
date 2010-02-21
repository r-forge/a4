\name{plotComb2Samples}
\alias{plotComb2Samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots the correlation in gene expression between two samples}
\description{
Plots the correlation in gene expression between two samples. Each dot represents
a gene, and the dots have a density-dependent coloring.
Genes with exceptional behavior can be highlighted by showing their gene symbol. 
}
\usage{
plotComb2Samples(object, x, y, trsholdX = NULL, trsholdY = NULL, probe2gene = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ExpressionSet object for the experiment}
  \item{x}{String containing the name of the first sample. This should be a 
  the name of a column in the \code{exprs} data of the \code{expressionSet} object.}
  \item{y}{String containing the name of the second sample. See \code{x}}
  \item{trsholdX}{Vector of two values specifying the X-axis thresholds within which
  genes should be highlighted by their gene symbol.}
  \item{trsholdY}{Vector of two values specifying the Y-axis thresholds within which
  genes should be highlighted by their gene symbol.}
  \item{probe2gene}{Boolean indicating whether the probeset should be translated to a gene symbol
    (used for the default title of the plot)}
  \item{\dots}{Possibility to add extra plot options. See \code{\link{par}}}
}

\references{ ~put references to the literature/web site here ~ }
\author{W. Talloen}

\seealso{\code{\link{plotCombMultSamples}}}
\examples{
if (require(ALL)){
  data(ALL, package = "ALL")
  ALL <- addGeneInfo(ALL)

  plotComb2Samples(ALL,"84004", "01003",
    trsholdX = c(10,12), trsholdY = c(4,6),
	xlab = "a B-cell", ylab = "a T-cell")
}
}
%% Add one or more standard keywords, see file 'KEYWORDS' in the
%% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
