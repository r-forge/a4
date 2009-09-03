\name{annotationTable}
\alias{annotationTable}
\title{Function to Create an annotationTable}
\description{
  This function takes data to be displayed as well as 
  data containing hyperlinks corresponding to displayed
  data and constructs an object of class annotationTable
}
\usage{
annotationTable(displayData, displayCols = NULL, hrefData = NULL)
}
\arguments{
  \item{displayData}{data frame containing data that is meant to be displayed
                     in a LaTeX table}
  \item{displayCols}{list of named character vectors (of length one) that function
                     as key-value pairs; the names (keys) correspond to columns for which
                     the hyperlinks should be generated whereas the strings (values)
                     indicate what kind of link should be produced based on the
                     corresponding column in the \code{displayData}. The values should be
                     one of \code{"EntrezId"} or \code{"GOId"}.}
  \item{hrefData}{data frame containing hyperlink information for the columns of
                  the same name in the \code{displayData} data frame}
}
\details{
  If \code{hrefData} is given, the \code{displayCols} are not taken into account.
  If no \code{hrefData} is given, the information in \code{displayCols} allows to 
  automatically create the \code{hrefData}.
}
\value{
  object of class 'annotationTable'
}
\author{Tobias Verbeke}
\examples{
  ## some dummy data
  dData <- data.frame(someSymbol = LETTERS[1:5], 
                      accessionNumber = 
                        c("X83928", "V00540", "U21090", "L38487", "M34057"))

  at <- annotationTable(displayData = dData, 
                        displayCols = list(accessionNumber = "EntrezId"))
}
\keyword{manip}
