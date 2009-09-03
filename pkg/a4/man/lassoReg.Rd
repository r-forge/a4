\name{lassoReg}
\alias{lassoReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multiple regression using the Lasso algorithm as implemented in the glmnet package}
\description{
Multiple regression using the Lasso algorithm as implemented in the glmnet package.
This is a theoretically nice approach to see which combination of genes predict best 
a continuous response. Empirical evidence that this actually works with high-dimensional
data is however scarce.
}
\usage{
lassoReg(object, covariate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object containing the expression measurements; currently the
    only method supported is one for ExpressionSet objects}
  \item{covariate}{character string indicating the column containing 
  the continuous covariate.}
}

\value{
  object of class \code{glmnet}
}
\references{ ~put references to the literature/web site here ~ }
\author{Willem Talloen}


\seealso{\code{\link{lassoClass}}}
\examples{
data(ALL, package = "ALL")
ALL <- addGeneInfo(ALL)
ALL$BTtype <- as.factor(substr(ALL$BT,0,1))

resultLasso <- lassoReg(object = ALL[1:100,], covariate = "age")
plot(resultLasso, label = TRUE,
		main = "Lasso coefficients in relation to degree of penalization.")
featResultLasso <- topTable(resultLasso, n = 15)
}


