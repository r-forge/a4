\name{volcanoPlot}
\alias{volcanoPlot}
\title{Draw a Volcano Plot}
\description{
  Generic function to draw a volcano plot. A volcano plot is a graph that allows to 
  simultaneously assess the P values (statistical significance) and log ratios
  (biological difference) of differential expression for the given genes.
}
\usage{
volcanoPlot(x, y, pointLabels, ...)
}
\arguments{
  \item{x}{either an object of class 'tTest', of class 'limma' or a numeric vector of 
   log ratios, i.e. the log of the fold change values; the names of the logRatio vector 
   will be used to display the names of the most interesting genes}
  \item{y}{should not be given if an object of class 'tTest' or 'limma' is passed as 
    argument 'x'; if 'x' is a numeric vector of log ratios, 'y' should be given and 
    should be a numeric vector of P-values indicating the statistical significance}
  \item{pointLabels}{Labels for points on the volcano plot that are interesting
    taking into account both the x and y dimensions; typically this is a
    vector of gene symbols; most methods can access the gene symbols directly from 
    the object passed as 'x' argument; the argument allows for custom labels if
    needed}
  \item{\dots}{further arguments to specific methods}
}
\value{
  The volcano plot is drawn to the current device.
}
\references{Goehlmann, H. and W. Talloen (2009). Gene Expression Studies Using Affymetrix
    Microarrays, Chapman \& Hall/CRC, pp. 148 - 153.
}
\author{Tobias Verbeke, based on code by Willem Talloen}
\examples{
data(ALL, package = "ALL")
ALL <- addGeneInfo(ALL)
ALL$BTtype <- as.factor(substr(ALL$BT,0,1))

tTestRes <- tTest(object = ALL,	groups = "BTtype", probe2gene = TRUE)
volcanoPlot(tTestRes)  
}
\keyword{dplot}
