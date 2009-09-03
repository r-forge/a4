\name{a4palette}
\alias{a4palette}
\title{Utility function that defines a color palette for use in a4}
\description{
  Utility function that defines a color palette for use in a4.
}
\usage{
a4palette(n)
}
\arguments{
  \item{n}{Number of color levels the palette should provide}
}
\details{
  For n = 1, \code{"darkblue"} is returned; for n = 2 
  \code{c("skyblue", "goldenrod")} is returned and for n > 2, 
  the output of \code{rainbow(n)} is returned.
}
\value{
  a character vector of colors
}
\author{Tobias Verbeke}
\seealso{\code{\link[grDevices]{rainbow}}}
\examples{
  a4palette(5)
}
\keyword{dplot}

