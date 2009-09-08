\name{probabilitiesPlot}
\alias{probabilitiesPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to plot the probabilities to belong to a certain class
 in binary classification problems.}
\description{
 Function to plot the probabilities to belong to a certain class
 in binary classification problems. These probabilities are often calculated 
 using a logistic regression model. The class membership of the samples is
  displayed using a colored strip (with legend below the plot).}
\usage{
probabilitiesPlot(proportions, classVar, sampleNames, plot = TRUE, barPlot = FALSE, layout = TRUE, main = NULL, sub = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{proportions}{A vector containing the calculated probabilities to belong to a certain class
 in binary classification problems. These probabilities are often calculated 
 using a logistic regression model.}
  \item{classVar}{A vector containing the class where the sample belongs to}
  \item{sampleNames}{A vector with the names of the samples}
  \item{plot}{logical.  If \code{FALSE}, nothing is plotted.}
  \item{barPlot}{Should a barplot be drawn (\code{TRUE}) or a scatterplot like
                 MCREstimate-type scores plot (the default, \code{FALSE}).}
  \item{layout}{boolean indicating whether \code{mcrPlot} should prespecify
    a layout for a single plot (default, \code{TRUE}) or whetherl the user
    takes care of the layout (\code{FALSE})}
  \item{main}{Main title for the scores plot; if not supplied, 'Scores Plot'
              is used as a default}
  \item{sub}{Subtitle for the scores plot; if not supplied, the classification
             technique and the chosen number of features are displayed}
  \item{\dots}{Additional graphical parameters to pass to the plot function}
}

\references{ ~put references to the literature/web site here ~ }
\author{Willem Talloen and Tobias Verbeke}

\seealso{\code{\link{logReg}}}
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
%% Add one or more standard keywords, see file 'KEYWORDS' in the
%% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
