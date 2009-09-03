\name{computeLogRatio}
\alias{computeLogRatio}
\title{ Summary statistics for gene expression}
\description{
Compute summary statistics per gene of expression data in a ExpressionSet object. 
}
\usage{
computeLogRatio(e, reference, within = NULL, across = NULL, nReplicatesVar = 3, ...)
}
\arguments{
  \item{e}{ An object of class ExpressionSet }
  \item{reference}{ A list with two items: var and level - See details }
  \item{within}{ Character vector - names of pData columns - See details }
  \item{across}{ Character vector - names of pData columns - See details  }
  \item{nReplicatesVar}{ Integer - Minimum number of replicates to compute variances }
  \item{\dots}{ \dots }
}
\details{
  Summary statistics (mean, variances and difference to reference or control) will be 
  computed on the 'exprs' slot of the ExpressionSet object. The parameters of the 
  computation are specified by the parameters 'reference', 'within' and 'across'.\\

  The design of the computations is such that the differences and pooled variances are calculated
  against the sample(s) that was(were) chosen as reference. The reference is specified by  
  the level of a certain variable in the phenoData slot (e.g.: column 'control' and level 'WT'
  of the phenoData slot or a boolean ('ref') variable with 0 or 1) -- the list object of 'var' and 
  'level' together determine the reference group. \\
  
  All groups determined by combining the \code{reference$var} and \code{across} variables will be 
  compared to the reference group. Two different approaches to obtain necessary computations: 
  \itemize{
    \item{}{Prepare a boolean variable that reflects only the reference group and specify all groupings 
       in the across arguments. E.g.: \code{reference=list(var = 'boolean', level = 1), 
       across = c('compound','dose')}}
    \item{}{Add an extra column to the phenoData slot that contains all combinations, with a specific one  
       for the reference group: for example, 
       \code{pData(e)['refvar'] <- paste(pData(e)['compound'], pData(e)['dose'],sep='.')} 
       so as to use \code{reference = list(var = 'refvar', level ='comp1.dose1')} as argument for reference.}
   } \\
   
  Sometimes computations need to be conducted within groups, and are thus nested. For example, 
  when comparing treament values of different cell lines, each will have gene expression values 
  for its own reference. The parameter 'within' allows to define such subgroups, for which 
  computations will be done separately and combined afterwards. Both parameters 'within' and 
  'across' can be a vector of column names, whose unique combinations will be used for groupings.    
}
\value{
  Returns an object of class ExpressionSet with pData inherited from the submitted ExpressionSet object, 
  supplemented by the computed statistics in the 'exprs' slot and info thereof in the 'phenoData' slot.
}
\author{Eric Lecoutre}
\seealso{\code{\link{plotLogRatio}}}
\examples{
data(ALL, package = "ALL")
ALL <- addGeneInfo(ALL)
ALL$BTtype <- as.factor(substr(ALL$BT,0,1))
ALL2 <- ALL[,ALL$BT != 'T1']  # omit subtype T1 as it only contains one sample
ALL2$BTtype <- as.factor(substr(ALL2$BT,0,1)) # create a vector with only T and B

# Test for differential expression between B and T cells
tTestResult <- tTest(ALL, "BTtype", probe2gene = FALSE)
topGenes <- rownames(tTestResult)[1:20]

# plot the log ratios versus subtype B of the top genes 
LogRatioALL <- computeLogRatio(ALL2, reference=list(var='BT',level='B'))
a <- plotLogRatio(e=LogRatioALL[topGenes,],openFile=FALSE, tooltipvalues=FALSE, device='X11',
		colorsColumnsBy=c('BTtype'), main = 'Top 20 genes most differentially between T- and B-cells',
		orderBy = list(rows = "hclust"),
		probe2gene = TRUE)
a <- plotLogRatio(e=LogRatioALL[topGenes,],openFile=TRUE, tooltipvalues=FALSE, device='pdf',
		colorsColumnsBy=c('BTtype'), main = 'Top 20 genes most differentially between T- and B-cells',
		orderBy = list(rows = "hclust", cols = "sex"),
		probe2gene = TRUE)
}
\keyword{manip}
\keyword{data}
\keyword{dplot}

