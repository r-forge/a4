\name{filterVarInt}
\alias{filterVarInt}
\title{Filter Features On Intensity and Variance}
\description{
  Function to filter on intensity and variance as
  typically used in gene expression studies
}
\usage{
filterVarInt(object, IntCutOff = log2(100), IntPropSamples = 0.25, VarCutOff = 0.5)
}
\arguments{
  \item{object}{ExpressionSet object}
  \item{IntCutOff}{cut-off value used for the intensity filter}
  \item{IntPropSamples}{proportion of samples used by the intensity
    filter; by default \code{IntPropSamples} si set to 0.25 is }
  \item{VarCutOff}{cut-off value used for the variance filter}
}
\details{
  The intensity filter implies that (by default) the intensity
  levels must be greater than log2(100) in at least 25 percent
  of the samples.
  
  The variance filter requires that the features have an interquartile
  range (IQR) greater than 0.5. Note that the IQR is quite insensitive 
  to outliers such that genes with outlying expression values in a few 
  samples are excluded as long as their overall variation is small.
}
\value{
  Object of class ExpressionSet containing only the features that
  pass the variance and intensity filter.
}
\references{
  Gentleman, R. et al. (2005). Bioinformatics and Computational Biology Solutions 
    using R and BioConductor, New York: Springer.
    
  Goehlmann, H. and W. Talloen (2009). Gene Expression Studies Using Affymetrix
    Microarrays, Chapman \& Hall/CRC, p. 128.
}
\author{Willem Talloen}
\seealso{\code{\link[genefilter]{pOverA}}, \code{\link[genefilter]{filterfun}}}
\examples{
  data(ALL, package = "ALL")
  fALL <- filterVarInt(ALL)
  fALL
}
\keyword{manip}