\name{volcanoPlot-methods}
\alias{volcanoPlot-methods}
\alias{volcanoPlot,tTest,missing,missing-method}
\alias{volcanoPlot,tTest,missing,character-method}
\alias{volcanoPlot,numeric,numeric,character-method}
\alias{volcanoPlot,numeric,numeric,missing-method}
\alias{volcanoPlot,limma,missing,missing-method}
\alias{volcanoPlot,limma,missing,character-method}

\title{Draw a Volcano Plot}
\description{
  This function draws a volcano plot, a graph that allows to simultaneously
  assess the statistical and biological significance of differential expression
  for the given genes.
}
\section{Methods}{
  \describe{

tTest,missing,missing

\item{x = "tTest", y = "missing", pointLabels = "missing"}{volcanoPlot
  for an object resulting from \code{tTest}}

tTest,missing,character

\item{x = "tTest", y = "missing", pointLabels = "character"}{volcanoPlot
  for an object resulting from \code{tTest}}

numeric,numeric,character

\item{x = "numeric", y = "numeric", pointLabels = "character"}{volcanoPlot
  for arbitrary numeric vectors containing log ratio values and p values
  respectively}

numeric,numeric,missing

\item{x = "numeric", y = "numeric", pointLabels = "missing"}{volcanoPlot
  for arbitrary numeric vectors containing log ratio values and p values
  respectively}

limma,missing,missing

\item{x = "limma", y = "missing", pointLabels = "missing"}{volcanoPlot
  for an object resulting from \code{limma2Groups}}
  
limma,missing,character

\item{x = "tTest", y = "missing", pointLabels = "missing"}{volcanoPlot
  for an object resulting from \code{limma2Groups}}

}
}
\arguments{
  \item{x}{either an object of class 'tTest', or a numeric vector of log ratios, i.e. the log of the fold change values;
    the names of the logRatio vector will be used to display the names of the most 
    interesting genes}
  \item{y}{should not be given if an object of class 'tTest' is passed as argument 'x';
    if 'x' is a numeric vector of log ratios, 'y' should be given and should be a 
    numeric vector of P-values indicating the statistical significance}
  \item{pointLabels}{Labels for points on the volcano plot that are interesting
    taking into account both the x and y dimensions; typically this is a
    vector of gene symbols; most methods can access the gene symbols directly from 
    the object passed as 'x' argument; the argument allows for custom labels if
    needed}
  \item{topPValues}{top n points that will be included in the points to label based
    on their low P Values}
  \item{topLogRatios}{top n points that will be included in the points to label based
    on their high absolute values of the log ratio}
  \item{smoothScatter}{use color saturation to indicate dots that are in densely
    populated regions of the graph; defaults to \code{TRUE}}
  \item{xlab}{label for the x axis (string)}
  \item{ylab}{label for the y axis (string)}
  \item{main}{main title for the graph (string)}
  \item{sub}{subtitle for the graph (string)}
}
\details{
  The set of genes for which labels are displayed is the \emph{union} of the set of
  genes that have lowest P-values (\code{topPValues}) and the set of genes
  that display the highest absolute values for the log ratios (\code{topLogRatios}).
}
\value{
  The volcano plot is drawn to the current device.
}
\author{Tobias Verbeke, based on code by Willem Talloen}
\keyword{dplot}
