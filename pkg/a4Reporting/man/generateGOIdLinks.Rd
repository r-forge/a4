\name{generateGOIdLinks}
\alias{generateGOIdLinks}
\title{Transform a GO ID into a hyperlink}
\description{
  Utility function to transform a GO ID into a hyperlink
  to the corresponding page on the gene ontology website
}
\usage{
generateGOIdLinks(x)
}
\arguments{
  \item{x}{vector of GO IDs}
}
\value{
  vector of hyperlinks
}
\author{Tobias Verbeke}
\note{ 
  Snippet taken from the \code{annaffy} package
}
\seealso{\code{\link{generateEntrezIdLinks}}}
\keyword{manip}
