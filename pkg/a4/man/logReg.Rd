\name{logReg}
\alias{logReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Logistic regression for predicting the probability to belong to a certain class
 in binary classification problems.}
\description{
Logistic regression for predicting the probability to belong to a certain class
 in binary classification problems.
}
\usage{
logReg(object, groups, probesetId = NULL, geneSymbol = NULL, main = NULL, probe2gene = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
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

\seealso{\code{\link{ROCcurve}}}
\examples{
data(ALL, package = "ALL")
ALL <- addGeneInfo(ALL)
ALL$BTtype <- as.factor(substr(ALL$BT,0,1))

logRegRes <- logReg(geneSymbol = "HLA-DPB1", object = ALL, groups = "BTtype")
ProbabilitiesPlot(proportions = logRegRes$fit, classVar = logRegRes$y,
		sampleNames = rownames(logRegRes), main = 'Probability of being a T-cell type ALL')
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
