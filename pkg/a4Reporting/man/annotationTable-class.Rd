\name{annotationTable-class}
\docType{class}
\alias{annotationTable-class}
\alias{show,annotationTable-method}

\title{Class "annotationTable"}
\description{Class to represent both displayed information and hyperlink information
  to prepare tabular output for LaTeX (with hyperlinks)}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("annotationTable", ...)}
or using the wrapper function \code{annotationTable}
}
\section{Slots}{
	 \describe{
    \item{\code{displayData}:}{Object of class \code{"data.frame"} containing the columns to be
      displayed in the table}
    \item{\code{displayCols}:}{Object of class \code{"list"} giving key-value pairs that allow
      to automatically generate the hyperlinks for the corresponding columns of the \code{displayData}}
    \item{\code{hrefData}:}{Object of class \code{"data.frame"} giving the hyperlink information for
      the corresponding columns of the \code{displayData}}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "annotationTable")}: print an annotationTable 
      (without displaying the hyperlink information)}
  }
}
\author{Tobias Verbeke}
\examples{
  showClass("annotationTable")
}
\keyword{classes}
