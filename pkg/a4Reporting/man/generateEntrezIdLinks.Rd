\name{generateEntrezIdLinks}
\alias{generateEntrezIdLinks}
\title{Transform an ENTREZ ID into a hyperlink}
\description{
  Utility function to transform an ENTREZ ID into a hyperlink on 
  the NCBI Entrez page for the given gene
}
\usage{
generateEntrezIdLinks(x)
}
\arguments{
  \item{x}{vector of Entrez IDs}
}
\value{
  vector of hyperlinks for the corresponding Entrez IDs
}
\author{Tobias Verbeke}
\note{ 
  Snippet taken from the \code{annaffy} package
}
\seealso{\code{\link{generateGOIdLinks}}}
\keyword{manip}
