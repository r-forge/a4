\name{profilesPlot}
\alias{profilesPlot}
\title{Plot expression profiles of multiple genes or probesets}
\description{
Plot expression profiles of multiple genes or probesets. Each line depicts a gene,
and the color legend can be used to identify the gene.
}
\usage{
profilesPlot(object, probesetIds, sampleIDs = TRUE, addLegend = TRUE, legendPos = "topleft", colvec = NULL, orderGroups = NULL, ...)
}
\arguments{
  \item{object}{ExpressionSet object for the experiment}
  \item{probesetIds}{The probeset ID. These should be stored in the \code{featureNames}
   of the \code{expressionSet} object.}
  \item{colvec}{Vector of colors to be used for the groups. If not specified, the default colors of
   \code{a4palette} are used.}
  \item{sampleIDs}{A boolean or a string to determine the labels on the x-axis. Setting it to FALSE
  results in no labels (interesting when the labels are unreadable due to large sample sizes).
  Setting it to a string will put the values of that particular \code{pData} column as labels.
  The string should be a name of a column in the \code{pData} of the \code{expressionSet} object."}
  \item{addLegend}{Boolean indicating whether a legend for the colors of the dots should be added.}
  \item{legendPos}{Specify where the legend should be placed. Typically either \code{topright},
   \code{bottomright}, \code{topleft} (the default) or \code{bottomleft}}
  \item{orderGroups}{String containing the name of the grouping variable to order the samples 
  in the x-axis accordingly. This should be a name of a column in the \code{pData} of
   the \code{expressionSet} object.}
  \item{\dots}{Possibility to add extra plot options. See \code{\link{par}}}
}

\author{W. Talloen}
\seealso{\code{\link{plot1gene}}, \code{\link{boxPlot}}}
\examples{
if (require(ALL)){
  data(ALL, package = "ALL")
  ALL <- addGeneInfo(ALL)
  ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
  
  myGeneSymbol <- c("LCK")	# a gene 
  probesetPos <- which(myGeneSymbol == featureData(ALL)$`Gene Symbol`)
  myProbesetIds <- featureNames(ALL)[probesetPos]
  
  # profilesPlot(object = ALL, probesetIds = myProbesetIds, 
  #     orderGroups = "BT", sampleIDs = "BT")
}
}
%% Add one or more standard keywords, see file 'KEYWORDS' in the
%% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
