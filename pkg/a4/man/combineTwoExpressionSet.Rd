\name{combineTwoExpressionSet}
\alias{combineTwoExpressionSet}
\title{ Combine two ExpressionSet objects }
\description{
  Merge two ExpressionSet objects, checking their attributes.
}
\usage{
combineTwoExpressionSet(x, y)
}
\arguments{
  \item{x}{ An object of class ExpressionSet }
  \item{y}{ An object of class ExpressionSet }
}
\details{
  exprs and pData are merged. Other data (such as MIAME or annotation) 
  are those of x.
}
\value{
  An object of class ExpressionSet
}
\author{Eric Lecoutre}
\seealso{  \code{\link[Biobase]{ExpressionSet}}}
\examples{
\dontrun{
# prepare and combine two ExpressionSet

data(data.H2009); data(phenoData.H2009)
data(data.SKOV3); data(phenoData.SKOV3)

eH2009 <- prepareExpressionSet(exprs = data.H2009, phenoData = phenoData.H2009, changeColumnsNames = TRUE)
eSKOV3  <- prepareExpressionSet(exprs = data.SKOV3, phenoData = phenoData.SKOV3, changeColumnsNames = TRUE)

newE <- combineTwoExpressionSet(eH2009,eSKOV3)
}
}
\keyword{data}

