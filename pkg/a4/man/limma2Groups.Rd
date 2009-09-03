\name{limma2Groups}
\alias{limma2Groups}
\title{Wrapper function for the comparison of two groups using limma}
\description{
  Wrapper function for the comparison of two groups using limma
}
\usage{
limma2Groups(object, group, probe2gene = TRUE)
}
\arguments{
  \item{object}{object of class ExpressionSet}
  \item{group}{string indicating the variable defining the two groups
               to be compared}
  \item{probe2gene}{logical; if \code{TRUE} Affymetrix probeset IDs are translated
    into gene symbols; if \code{FALSE} no such translation is done}
}
\value{
  S4 object of class 'limma' with the following two components:
  \item{MArrayLM}{S4 object of class MArrayLM as returned by the limma
    function of the limma package}
  \item{geneSymbols}{character vector of gene symbols; this slot is only
    populated if \code{probe2gene=TRUE} (and if the ExpressionSet object
    is appropriately annotated by \code{addGeneInfo} for gene symbols to
    be extracted)}
}
% TODO add references to book
% \references{}
\author{Tobias Verbeke and Willem Talloen}
\note{
  A 'topTable' method is defined for 'limma' objects.
}
% \seealso{}
% \examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

