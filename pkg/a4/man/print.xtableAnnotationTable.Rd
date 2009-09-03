% TODO add examples section
\name{print.xtableAnnotationTable}
\alias{print.xtableAnnotationTable}
\title{Print method for 'xtableAnnotationTable' objects}
\description{
  Print or export 'xtableAnnotationTable' objects
}
\usage{
print.xtableAnnotationTable(x, ...)
}
\arguments{
  \item{x}{Object of class 'xtableAnnotationTable'}
  \item{\dots}{Further arguments passed to \code{print.xtable}}
}
\details{
  Wrapper to be able to use a specific \code{sanitize.text} function
  in the \code{print.xtable} call
}
\author{Tobias Verbeke}
\seealso{\code{\link[xtable]{print.xtable}}}
\keyword{manip}
