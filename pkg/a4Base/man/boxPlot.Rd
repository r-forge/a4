\name{boxPlot}
\alias{boxPlot}
\title{Create a boxplot for a given gene.}
\description{
  Create a boxplot for a given gene. The boxplot displays the expression values (y-axis)
  by groupss (x-axis). The raw data are superimposed as dots, jittered for readability of the plot.
  Optionally, the dots can be colored by another variable.
}
\usage{
boxPlot(probesetId = NULL, geneSymbol = NULL, object, groups, main = NULL, colvec = NULL, colgroups = NULL, probe2gene = TRUE, addLegend = TRUE, legendPos = "topleft", ...)
}
\arguments{
  \item{probesetId}{The probeset ID. These should be stored in the \code{featureNames}
   of the \code{expressionSet} object.}
  \item{geneSymbol}{The gene symbol. These should be stored in the column \code{`Gene Symbol`}
   in the \code{featureData} of the \code{expressionSet} object.}
  \item{object}{ExpressionSet object for the experiment}
  \item{groups}{String containing the name of the grouping variable. This should be a 
  the name of a column in the \code{pData} of the \code{expressionSet} object.}
  \item{main}{Main title on top of the graph}
  \item{colvec}{Vector of colors to be used for the groups. If not specified, the default colors of
   \code{a4palette} are used.}
  \item{colgroups}{String containing the name of the variable to color the superimposed dots.
   This should be a the name of a column in the \code{pData} of the \code{expressionSet} object.}
  \item{probe2gene}{Boolean indicating whether the probeset should be translated to a gene symbol
    (used for the default title of the plot)}
  \item{addLegend}{Boolean indicating whether a legend for the colors of the dots should be added.}
  \item{legendPos}{Specify where the legend should be placed. Typically either \code{topright},
   \code{bottomright}, \code{topleft} (the default) or \code{bottomleft}}
  \item{\dots}{Possibility to add extra plot options. See \code{\link{par}}}
}


\author{Willem Talloen}

\seealso{\code{\link{plot1gene}}}
\examples{
# simulated data set
esSim <- simulateData()
boxPlot(probesetId = 'Gene.1', object = esSim, groups = 'type', addLegend = FALSE)

# ALL data set
if (require(ALL)){
  data(ALL, package = "ALL")
  ALL <- addGeneInfo(ALL)
  ALL$BTtype <- as.factor(substr(ALL$BT,0,1))

  boxPlot(geneSymbol = 'HLA-DPB1', object = ALL, boxwex = 0.3,
		  groups = 'BTtype', colgroups = 'BT', legendPos='topright')
}
}
