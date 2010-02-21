\name{volcanoplotter}
\alias{volcanoplotter}
\title{Workhorse function for the different volcanoPlot methods}
\description{
  Workhorse function for the different volcanoPlot methods. A volcano plot 
  is a graph that allows to simultaneously assess the P values (statistical 
  significance) and log ratios (biological difference) of differential 
  expression for the given genes.
}
\usage{
volcanoplotter(logRatio, pValue, pointLabels, topPValues = 10, topLogRatios = 10, logTransformP = TRUE, smoothScatter = TRUE, xlab = NULL, ylab = NULL, main = NULL, sub = NULL, newpage = TRUE)
}
\arguments{
  \item{logRatio}{numeric vector of log ratios}
  \item{pValue}{numeric vector of P values}
  \item{pointLabels}{Labels for points on the volcano plot that are interesting
    taking into account both the x and y dimensions; typically this is a
    vector of gene symbols; most methods can access the gene symbols directly from 
    the object passed as 'x' argument; the argument allows for custom labels if
    needed}
  \item{topPValues}{top n points that will be included in the points to label based
    on their low P Values}
  \item{topLogRatios}{top n points that will be included in the points to label based
    on their high absolute values of the log ratio}
  \item{logTransformP}{if \code{TRUE} (default) -log10(pValue) is used for the plot 
    instead of the raw P values}
  \item{smoothScatter}{use color saturation to indicate dots that are in densely
    populated regions of the graph; defaults to \code{TRUE}}
  \item{xlab}{label for the x axis (string)}
  \item{ylab}{label for the y axis (string)}
  \item{main}{main title for the graph (string)}
  \item{sub}{subtitle for the graph (string)}
  \item{newpage}{should the graph be drawn to a new grid page? Defaults to
    \code{FALSE}. This argument is useful for including several volcano plots 
    in one layout.}
}
\value{
  a volcanoplot is drawn to the current device
}
\author{Tobias Verbeke}
\seealso{\code{\link{volcanoPlot-methods}}}
\keyword{dplot}
