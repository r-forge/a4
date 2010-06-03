\name{xtable-methods}
\docType{methods}
\alias{xtable-methods}
\alias{xtable,annotationTable,missing,missing,missing,missing,missing-method}
\alias{xtable,annotationTable,ANY,ANY,ANY,ANY,ANY-method}
\alias{xtable,annotationTable,ANY,ANY,ANY,numeric,ANY-method}
\alias{xtable.topTableGlmnet}
\alias{xtable.pamClassConfusionTable}
\alias{xtable.topTablePam}
\title{Methods for Function xtable in Package `annotationTable'}
\description{
   \code{xtable} methods for several a4 objects, such as \code{annotationTable} objects, 
   \code{topTable} objects etc.
}
\section{Methods}{
\describe{

\item{x = "annotationTable", caption = "missing", label = "missing", align = "missing", digits = "missing", display = "missing"}{generates a LaTeX representation for the given annotationTable}
\item{x = "annotationTable", caption = "ANY", label = "ANY", align = "ANY", digits = "ANY", display = "ANY"}{generates a LaTeX representation for the given annotationTable}
\item{x = "annotationTable", caption = "ANY", label = "ANY", align = "ANY", digits = "numeric", display = "ANY"}{generates a LaTeX representation for the given annotationTable}

}}
\examples{
  ## some dummy data
  dData <- data.frame(someSymbol = LETTERS[1:5], 
                      accessionNumber = 
                        c("X83928", "V00540", "U21090", "L38487", "M34057"))

  at <- annotationTable(displayData = dData, 
                        displayCols = list(accessionNumber = "EntrezId"))
  xat <- xtable(at)
  print(xat, include.rownames = FALSE)
}
\keyword{methods}
\keyword{manip}
