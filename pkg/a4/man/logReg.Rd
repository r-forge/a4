\name{logReg}
\alias{logReg}
\title{Logistic regression for predicting the probability to belong to a certain class
 in binary classification problems.}
\description{
Logistic regression for predicting the probability to belong to a certain class
 in binary classification problems.
}
\usage{
logReg(object, groups, probesetId = NULL, geneSymbol = NULL, main = NULL, probe2gene = TRUE, ...)
}
\arguments{
  \item{object}{ExpressionSet object for the experiment}
  \item{groups}{String containing the name of the grouping variable. This should be a 
  the name of a column in the \code{pData} of the \code{expressionSet} object.}
  \item{probesetId}{The probeset ID. These should be stored in the \code{featureNames}
   of the \code{expressionSet} object.}
  \item{geneSymbol}{The gene symbol. These should be stored in the column \code{`Gene Symbol`}
   in the \code{featureData} of the \code{expressionSet} object.}
  \item{main}{Main title on top of the graph}
  \item{probe2gene}{Boolean indicating whether the probeset should be translated to a gene symbol
    (used for the default title of the plot)}
  \item{\dots}{Possibility to add extra plot options. See \code{\link{par}}}
}
\details{
It will always estimate probability scores to belong to the second level
of the factor variable. If a probability score to other level is preferred,
then you need to change the order of the levels of the factor.
}
\value{
A data.frame object with three columns and rownames
  \item{rownames }{The 'sampleNames' of the expressionSet object}
  \item{x }{The expression values for the specified gene for all samples}
  \item{y }{The labels of the samples}
  \item{fit }{The fitted probability score to belong to one of the two classes.}
}
\references{ ~put references to the literature/web site here ~ }
\author{Willem Talloen}

\seealso{\code{\link{ROCcurve}},\code{\link{probabilitiesPlot}}}
\examples{
data(ALL, package = "ALL")
ALL <- addGeneInfo(ALL)
ALL$BTtype <- as.factor(substr(ALL$BT,0,1))

logRegRes <- logReg(geneSymbol = "HLA-DPB1", object = ALL, groups = "BTtype")

# scoresplot
probabilitiesPlot(proportions = logRegRes$fit, classVar = logRegRes$y,
		sampleNames = rownames(logRegRes), main = 'Probability of being a T-cell type ALL')
# barplot
probabilitiesPlot(proportions = logRegRes$fit, classVar = logRegRes$y, barPlot=TRUE,
		sampleNames = rownames(logRegRes), main = 'Probability of being a T-cell type ALL')
}
