\name{createExpressionSet}
\alias{createExpressionSet}
\title{ combine gene expression and phenotype data onto a ExpressionSet object}
\description{
  Basically a wrapper for \code{new('ExpressionSet',...)}, this function gathers gene
  expression and phenotype data, after having checked their compatibility.
}
\usage{
createExpressionSet(exprs = new("matrix"), phenoData = new("AnnotatedDataFrame"), varMetadata = NULL, dimLabels = c("rowNames", "colNames"), featureData = NULL, experimentData = new("MIAME"), annotation = character(0), changeColumnsNames = TRUE, ...)
}
\arguments{
  \item{exprs}{ gene expression matrix }
  \item{phenoData}{ phenotype data associated with exprs columns, as a matrix or data.frame }
  \item{varMetadata}{ optional metadata on phenotype data}
  \item{dimLabels}{ see 'ExpressionSet'}
  \item{featureData}{ see 'ExpressionSet' }
  \item{experimentData}{ see 'ExpressionSet' }
  \item{annotation}{ see 'ExpressionSet' }
  \item{changeColumnsNames}{ Change exprs columns names -- see details }
  \item{\dots}{ \code{\dots} }
}
\details{
   If \code{changeColumnsNames} is \code{TRUE}, then the procedure is the following: first one checks if phenoData contains a column named 'colNames'. If so, content will be used to rename exprs colums. On the other case, one uses combinations of phenoData columns to create new names. In any case, old columns names
   are stored within a column named 'oldcolnames' in the pData.
}
\value{
  An object of class ExpressionSet
}
\author{Eric Lecoutre}
\seealso{ \code{\link[Biobase]{ExpressionSet}}}
\examples{
# simulate expression data of 10 features (genes) measured in 4 samples
x <- matrix(rnorm(40), ncol = 4)
colnames(x) <- paste("sample", 1:4, sep = "_")
rownames(x) <- paste("feature", 1:10, sep = "_")

# simulate a phenodata with two variables
ToBePheno <- data.frame(Gender = rep(c('Male','Female'), 2), 
		Treatment = rep(c('Trt','Control'), each=2))
rownames(ToBePheno) <- paste("sample", 1:4, sep = "_")

eset <- createExpressionSet(exprs = x, phenoData = ToBePheno)
}

\keyword{data}
