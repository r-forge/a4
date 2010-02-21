\name{plotLogRatio}
\alias{plotLogRatio}
\title{ Plot a summary gene expression graph }
\description{
Plot ratios of expression values observed in a treatment versus those of a reference. 
First the ratios and variances are computated on the gene expression data.
}
\usage{
plotLogRatio(e, reference, within = NULL, across = NULL, nReplicatesVar = 3,
 filename = "Rplots", device = "svg", orderBy = list(rows = "hclust", cols = NULL),
  colorsColumns = NULL, colorsColumnsBy = NULL, colorsColumnsByPalette = c("#1B9E77", "#D95F02",
   "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"), colorsUseMeanQuantiles = FALSE,
    colorsMeanQuantilesPalette = c("orange", "red", "darkred"), colorsBarsMatrix = NULL,
     colorsGenesNames = c("black"), main = paste("log2 ratio's"), shortvarnames = NULL,
      longvarnames = NULL, gene.length = 50, gene.fontsize = 6, main.fontsize = 9, 
      columnhead.fontsize = 8, mx = 1.5,
       exp.width = 1.8, exp.height = 0.2, log2l.show = TRUE, log4l.show = FALSE, quantiles.show = FALSE,
        quantiles.compute = c(0.9), error.show = TRUE, view.psid = FALSE,
         errorLabel = "Error bars show the pooled standard deviation", 
          closeX11 = FALSE, openFile = FALSE, tooltipvalues = FALSE, probe2gene = TRUE, ...)
}
\arguments{
  \item{e}{ ExpressionSet object to use }
  \item{reference}{ List with components 'var' and 'level ' -- see computeDEStats help }
  \item{within}{ Vector of characters for pData column -- see computeDEStats help }
  \item{across}{ Vector of characters for pData column -- see computeDEStats help }
  \item{nReplicatesVar}{ Minimum number of replicates to compute variances and pooled 
    standard errors -- see computeDEStats help }
  \item{filename}{ Name of the filename to use. No need to specify extension which 
    will be added according to device.}
  \item{device}{ One of 'pdf', 'X11', 'png', 'svg'. For svg device, one X11 device 
    is also opened.}
  \item{orderBy}{ See details }
  \item{colorsColumns}{A vector of colors to be used for plotting columns; default value 
    is NULL which ends up with red -- see Colors section}
  \item{colorsColumnsBy}{ A vector of pData columns which combinations specify different colors to be used -- see Colors section}
  \item{colorsColumnsByPalette}{If colorsColumns is NULL, vector of colors to be used for coloring columns potentially 
  splitted by colorsColumnsBy}
  \item{colorsUseMeanQuantiles}{Boolean to indicate if the quantile groups computed 
    on averages over all treatments should be used for coloring -- see Colors section}
  \item{colorsMeanQuantilesPalette}{if colorsUseMeanQuantiles is TRUE, these colors 
    will be used for the different groups -- see Colors section}
  \item{colorsBarsMatrix}{Matrix of colors to be used for each individual bar; 
    colors are provided for genes in data order and thus are possibly reordered 
    according to orderBy -- see Colors section}
  \item{colorsGenesNames}{Vector of colors to be used for gene names; will be recycled 
    if necessary; colors are provided for genes in data order and thus are possibly 
    reordered according to orderBy}
  \item{main}{Main title}
  \item{shortvarnames}{Vector or pData column to be used to display in graph columns. If NULL, those names 
    will be used from the coded names added to pData during computations (list of columns values pasted with a dot). Warning: shortvarnames must be defined in the order columns are present in the ExpressionSet object so that they will be reordered if one asks to order columns.}
  \item{longvarnames}{ pData column to be used in SVG tooltip title. If NULL, shortvarnames will be used.Same warning than shortvarnames about ordering. }
  \item{gene.length}{ Maximum number of characters that will be printed of the gene names}
  \item{gene.fontsize}{ Font size for the gene names , default = 6}
  \item{main.fontsize}{ Font size for the main, default = 9}
  \item{columnhead.fontsize}{ Font size for the column headers, default = 8}
 
  \item{mx}{ Expansion factor for the width of the bars that represent the expression ratios}
  \item{exp.width}{ Expansion factor for global graph width, and the space between the plotted columns }
  \item{exp.height}{ Expansion factor for global graph height, and the space between the plotted rows }
  \item{log2l.show}{ A logical value. If 'TRUE', the line for log2 values on each column (when max(data) > 2) is drawn}
  \item{log4l.show}{ A logical value. If 'TRUE', the line for log4 values on each column (when max(data) > 4) is drawn}
  \item{quantiles.show}{ A logical value. If 'TRUE', a line is drawn for quantiles computed separately on each columns }
  \item{quantiles.compute}{ A logical value. If 'TRUE', the vector quantiles will be computed and displayed
    provided that \code{quantile.show} is \code{TRUE}}
  \item{error.show}{ A logical value. If 'TRUE', errors bars are displayed on the graph 
    (only for those columns for which they are available) }
  \item{view.psid}{ A logical value. If 'TRUE', the genes psid is displayed on the gene names }
  \item{errorLabel}{ A character vector describing the error bars, printed at the bottom of the figure}
  \item{closeX11}{If \code{device} is SVG, do we close the required X11 device at the end?}
  \item{openFile}{ A logical value. If 'TRUE', the produced output file is opened }
  \item{tooltipvalues}{ If device is SVG, one can choose to display each bar separately, with data values as tooltips. 
    Note however that each bar will be considered as a distinct object instead of a column, which will takes much 
    more time to create the graph and produces a much bigger SVG file}
  \item{probe2gene}{Boolean indicating whether the probeset should be translated to a gene symbol
    (used for the default title of the plot)}
  \item{\dots}{ \code{\dots} }
}
\section{Ordering}{
orderBy: A list with two components, rows and cols, each one possibly being NULL (no ordering 
on the specific dimension). Ordering on cols can be done according to (a) pData column(s) 
(for example: \code{c('cellline','compound','dose'}. Ordering on rows can be done using of the 
following values:
 \itemize{
  \item{NULL}{no reordering on rows}
  \item{numeric vector}{use the vector values to sort rows}
  \item{alpha}{use genes names alphabetice order}
  \item{effect}{try to assess global gene expression level by taking sum(abs(values)) on specified exprs columns)}
  \item{hclust}{use the ordering returned by \code{hclust} invoked on specified exprs columns}
 }
 }
\section{Colors}{
  The management of colors is very flexible but is a little bit tricky, as a variety of parameters 
  are available to the user. Basically, combinations of arguments allow to set colors for columns headers (text), 
  columns as a whole (different colors for the different columns) or for each of the inividual horizontal bars.
  By default, everything is red. There are four main different arguments that can be used and that are 
  applied in a consecutive order. Each one may override a previous argument value. Below is a list of 
  arguments and their consecutive actions:
 \itemize{
  \item{\code{colorsColumns}}{ The first way to assign colors is to provide a vector of colors that will 
    be used for each column (headers and its horizontal bars). This vector is recycled so that providing one unique 
    value will color all columns, whereas providing a vector of length 2 will alternate columnns colors.}
  \item{\code{colorsColumnsBy}}{To be used when the experiment involves groupings for pData, for example dose, 
    cellline or treatment. In order to see the effects of such variables, one can color columns using 
    combinations of those. The argument is a vector of pData columns such as \code{c('cellline','dose')}. 
    Unique combinations will be computed and a color will be assigned for each group of columns. 
    The vector that is provided with the argument \code{colorsColumnsByPalette} is used to assign colors. 
    If the argument \code{colorColumnsBy} is not \code{NULL} then it overrides the previous argument \code{colorsColumns}.}
  \item{\code{colorsUseMeanQuantiles}}{ A logical value. The default plotGeneDE displays for each gene the expression value difference 
    between treatment and reference, but does not reveal any information about the expression levels in these conditions. 
    Parameter \code{colorsUseMeanQuantiles} allows to color the horizontal bars according to expression level that 
    is derived from quantiles computed on averages of the complete ExpressionSet object. 
    As it involves the expression data of all probesets, computations must be done 
    before subsetting the ExpressionSet object and the plotGeneDEting. The function \code{\link{addQuantilesColors}} 
    computes quantiles and corresponding mean expression level intervals. If \code{colorsUseMeanQuantiles} 'TRUE', 
    previous coloring parameters are overriden. The parameter \code{colorsMeanQuantilesPalette} is used to assign 
    colors for average-quantiles-groups. Note that columns headers are still given by previous arguments.}
  \item{\code{colorsBarsMatrix}}{The most flexible way to assign colors as the matrix will be used to color each bar 
    of the plot individually. A check is done to ensure that the number of rows and columns are not less than the number of 
    probesets and columns. If not \code{NULL}, this parameter overrides the previous ones.}
 }
 }  
\value{
  The ExpressionSet object with the computated variables is returned.
}
\author{Hinrich Goehlmann and Eric Lecoutre}
\seealso{ \code{\link{computeLogRatio}},\code{\link{addQuantilesColors}}}
\examples{
if (require(ALL)){
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
      orderBy = list(rows = "hclust"), probe2gene = TRUE)
}
\dontrun{		
  a <- plotLogRatio(e=LogRatioALL[topGenes,],openFile=TRUE, tooltipvalues=FALSE, device='pdf',
      colorsColumnsBy=c('BTtype'), main = 'Top 20 genes most differentially between T- and B-cells',
      orderBy = list(rows = "hclust", cols = "sex"), probe2gene = TRUE)
}
}
\keyword{hplot}
