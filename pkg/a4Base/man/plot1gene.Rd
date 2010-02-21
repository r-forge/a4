\name{plot1gene}
\alias{plot1gene}
\title{Create a Profile Plot for a given Gene}
\description{
  Create a profile plot for a given gene. A profile plot displays the expression values (y-axis)
  by samples (x-axis), sorted by group. This is a useful working graph as samples can be
  directly identified. For presentation purposes, a \code{boxPlot} can also be considered. with jittered for readability of the plot.
}
\usage{
plot1gene(probesetId = NULL, geneSymbol = NULL, object, groups, main = NULL, colvec = NULL,
	colgroups = NULL, probe2gene = TRUE, sampleIDs = TRUE, 
	addLegend = TRUE, legendPos = "topleft", ...)
}
\arguments{
  \item{probesetId}{The probeset ID. These should be stored in the \code{featureNames}
   of the \code{expressionSet} object.}
  \item{geneSymbol}{The gene symbol. These should be stored in the column \code{`Gene Symbol`}
   in the \code{featureData} of the \code{expressionSet} object.}
  \item{object}{ExpressionSet object for the experiment}
  \item{groups}{String containing the name of the grouping variable. This should be a 
  name of a column in the \code{pData} of the \code{expressionSet} object.}
  \item{colgroups}{String containing the name of the variable to color the superimposed dots.
   This should be a the name of a column in the \code{pData} of the \code{expressionSet} object.}
  \item{main}{Main title on top of the graph}
  \item{colvec}{Vector of colors to be used for the groups. If not specified, the default colors of
   \code{a4palette} are used.}
  \item{probe2gene}{Boolean indicating whether the probeset should be translated to a gene symbol
    (used for the default title of the plot)}
  \item{sampleIDs}{A boolean or a string to determine the labels on the x-axis. Setting it to FALSE
  results in no labels (interesting when the labels are unreadable due to large sample sizes).
  Setting it to a string will put the values of that particular \code{pData} column as labels.
  The string should be a name of a column in the \code{pData} of the \code{expressionSet} object."}
  \item{addLegend}{Boolean indicating whether a legend for the colors of the dots should be added.}
  \item{legendPos}{Specify where the legend should be placed. Typically either \code{topright},
   \code{bottomright}, \code{topleft} (the default) or \code{bottomleft}}
  \item{\dots}{Possibility to add extra plot options. See \code{\link{par}}}
}
\value{
  If a \code{geneSymbol} is given that has more than one probeSet,
  the plots for only the first probeSet is displayed.
  A character vector of corresponding probeset IDs is returned invisibly,
  so that one can check the profiles of the other related probeset IDs with
  an extra \code{plot1gene} statement
  
  If a \code{probesetId} is given, one single profile plot for the probeset is 
  displayed.
}
\references{ ~put references to the literature/web site here ~ }
\author{S. Osselaer, W. Talloen, T. Verbeke}
\seealso{\code{\link{plotCombination2genes}}, \code{\link{boxPlot}}}

\examples{
if (require(ALL)){
  data(ALL, package = "ALL")
  ALL <- addGeneInfo(ALL)

  # one variable (specified by groups)
  aa <- plot1gene(geneSymbol = 'HLA-DPB1', object = ALL, groups = "BT",
	    addLegend = TRUE, legendPos = 'topright')
  aa
  
  # two variables (specified by groups and colGroups)
  ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
  plot1gene(probeset = '1636_g_at', object = ALL, groups = 'BT',
      colgroups = 'mol.biol', legendPos='topright', sampleIDs = 'BT')    
}
}
\keyword{dplot}
