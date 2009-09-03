\name{plotCombMultSamples}
\alias{plotCombMultSamples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots the correlation in gene expression between more than 2 samples}
\description{
Plots a correlation matrix in gene expression between two samples. Each dot represents
a gene, and the dots have a density-dependent coloring.
}
\usage{
plotCombMultSamples(exprsMatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exprsMatrix}{ExpressionSet object to plot. For larger datasets,
   this will typically be a subset of the data.}
}

\references{ ~put references to the literature/web site here ~ }
\author{Willem Talloen}

\seealso{\code{\link{plotComb2Samples}}}
\examples{
data(ALL, package = "ALL")
ALL <- addGeneInfo(ALL)

plotCombMultSamples(exprs(ALL)[,c("84004", "11002", "01003")])
}
%% Add one or more standard keywords, see file 'KEYWORDS' in the
%% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
