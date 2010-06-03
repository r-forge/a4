\name{heatmap.expressionSet}
\alias{heatmap.expressionSet}

\title{
Image plot of an expressionSet
}

\description{
	Grid version of heatmap function adapted to expressionSet objects with some specific 
	requirements such as the possibility to display subgroups, define colors, 
	adapt text graphical parameters (sizes...).
	
	The function also suggests a size appropriate for a device to generate a complete plot with all elements. 
}

\usage{
heatmap.expressionSet(eset, col.groups = pData(phenoData(eset))[, "subGroup"], col.orderBy = order(pData(phenoData(eset))[, "subGroup"]), col.groups.sep.width = unit(8, "points"), col.labels = sampleNames(eset), col.labels.sep.width = unit(10, "points"), col.labels.gpar = gpar(cex = 1), col.labels.max.nchar = 20, colors.pergroup = FALSE, colors.groups = NULL, colors.groups.min = rgb(1, 1, 1), colors.max = rgb(1, 0, 0), colors.min = rgb(1, 1, 1), colors.nbreaks = 128, colors.palette = NULL, cell.gpar = gpar(lty = 0), row.groups.sep.height = unit(15, "points"), row.labels.sep.height = unit(10, "points"), row.col.groups.display = ifelse(length(unique(col.groups)) > 1, TRUE, FALSE), row.col.groups.display.height = unit(6, "points"), row.labels.gpar = gpar(cex = 1, col = "black"), row.labels.max.nchar = 45, row.labels = list("SYMBOL", "GENENAME"), row.labels.sep = " - ", row.groups = rep(1, nrow(exprs(eset))), row.order = "none", row.groups.hclust = FALSE, row.groups.hclust.n = 4, distfun = dist, hclustfun = function(d) {
    hclust(d, method = "ward")
}, values.min = 0, values.max = 16, title.gpar = gpar(cex = 1.4), title.main = "This is the title  possibly being very long - it will be splited on several lines or even displayed with dots at the end -- see there (does it work? addendum)", title.just = c("right", "top"), title.maxlines = 4, title.cutpoint = 40, subtitle.gpar = gpar(cex = 1), subtitle.main = "This is subtitle", subtitle.maxlines = 4, subtitle.just = title.just, subtitle.cutpoint = 40, margin.top = unit(2, "lines"), margin.left = unit(2, "lines"), margin.right = unit(2, "lines"), margin.bottom = unit(2, "lines"), legend.display = TRUE, legend.range = "full", legend.data.display = ifelse(legend.range == "full", TRUE, FALSE), legend.gpar = gpar(cex = 1), legend.width = unit(250, "points"), legend.height = unit(40, "points"), ...)
}


\arguments{
  \item{eset}{
		expressionSet object
}
  \item{col.groups}{
		Vector specifying sub-groups for individual. Sub-groups are treated separately and can thus on plot have different colors.
}
  \item{col.orderBy}{
		Vector specifying ordering for individual. In case there are sub-groups, individual must first be ordered by sub-groups, but an additional variable gives a way to sort individual within sub-groups. 
}
  \item{col.groups.sep.width}{
		Object of class unit (grid package). Width used to visually separate sub-groups of individuals. This can be unit(0,"points") for example for no separation.
}
  \item{col.labels}{
 	Character vector for columns labels (individuals), by default taken from phenoData.	
}
  \item{col.labels.sep.width}{
		Object of class unit (grid package). Space between image matrix zone and columns labels. 
}
  \item{col.labels.gpar}{
		Object of class gpar (grid package). Parameters to be used for labels (cex,...).
}
  \item{col.labels.max.nchar}{
		Integer. Number of maximum characters to be used for labels truncation.
}
  \item{colors.pergroup}{
		Boolean. If TRUE, separate colors are used to color image matrix. Colors defined for groups are used.
}
  \item{colors.groups}{
		Vector. Colors to be used for each group of individual. If NULL (default), colors are taken from column "sampleColor" of expressionSet phenodata.
}
  \item{colors.groups.min}{
		Character vector of length 1 corresponding to a valid color. If colors.groups are provided, a shading if done between color.group and this color (default: white). 
}
  \item{colors.max}{
		 Character vector of length 1 corresponding to a valid color. See colors details.
}
  \item{colors.min}{
		Character vector of length 1 corresponding to a valid color. See colors details.
}
  \item{colors.nbreaks}{
		Integer. Number of cutpoints used to split the color palette/shading. 
}
  \item{colors.palette}{
		Character vector of valid color names.
}
  \item{cell.gpar}{
		Object of class gpar (grid package). Parameters used to format cells, for example to add border (gpar(lty=1)).
}
  \item{row.groups.sep.height}{
		Object of class unit (grid package). Height between rows sub-groups.
}
  \item{row.labels.sep.height}{
		Object of class unit (grid package). Height between image plot zone and rows labels
}

  \item{row.col.groups.display}{
		Boolean. Display or not colored band for subgroups of individuals.
}
  \item{row.col.groups.display.height}{
		Object of class unit (grid package). If row.col.groups.display is TRUE then height used for the displayed band. 
}
  \item{row.labels.gpar}{
		Object of class gpar (grid package). Parameters to be used for labels (cex,...).
}
  \item{row.labels.max.nchar}{
		Integer. Number of maximum characters to be used for labels truncation.
}
  \item{row.labels}{
		Character vector or list. If vector, direct labels to be used. 
		If list, elements of the list will be taken from featureData and collapsed using row.labels.sep.%%  Can be a list    ~~Describe \code{row.labels} here~~
}
  \item{row.labels.sep}{
	In case labels are taken from featureData (list for row.labels), separator used to paste the provided columns.
}
  \item{row.groups}{
		Boolean specifying whether rows are split into sub-groups. %%     ~~Describe \code{row.groups} here~~
}
  \item{row.order}{
		Either a vector of indices to be used to reorder features (rows) or "none" or "hclust" to use clustering.
}
  \item{row.groups.hclust}{
	Boolean. If row.order equals "hclust", one can ask to split features into sub-groups based on a cut of the clustering dendogram.%%     ~~Describe \code{row.groups.hclust} here~~
}
  \item{row.groups.hclust.n}{
		Integer. If row.order equals "hclust" and row.groups.hclust is TRUE, number of sub-groups. 
}
  \item{distfun}{
		Function. For row.order equals "hclust", metric function.
}
  \item{hclustfun}{
		Function. For row.order equals "hclust", clustering function.
}
  \item{values.min}{
	Minimum value for the data range. Values that are inferior are assigned to that value. That ensures a maximal cutpoint for the coloring scale.
}
  \item{values.max}{
	Maximum value for the data range. Values that are superior are assigned to that value. That ensures a maximal cutpoint for the coloring scale.
}
  \item{title.gpar}{
		Object of class gpar (grid package). Parameters to be used for the main title (cex,...).
}
  \item{title.main}{
	Character vector. Main title to be displayed.
}
  \item{title.just}{
		Title justification, one of "center","left","right" (first letter of the word can also be used).  
}
  \item{title.maxlines}{
		Maximum number of lines for the title split.
}
  \item{title.cutpoint}{
		Integer. Maximum number of characters a line must have. Title is split into lines according to that cutpoint.}

  \item{subtitle.gpar}{
		Object of class gpar (grid package). Parameters to be used for the subtitle (cex, col,...).
}
  \item{subtitle.main}{
		Character vector. Subtitle. The subtitle will be split into lines following same rules as used for main title.
}
  \item{subtitle.maxlines}{
		Maximum number of lines for the subtitle split.
}
  \item{subtitle.just}{
		Subtitle justification, one of "center","left","right" (first letter of the word can also be used).  
}
  \item{subtitle.cutpoint}{
		Integer. Maximum number of characters a line must have. Subtitle is split into lines according to that cutpoint.%%Number of characters 
}
  \item{margin.top}{
  Object of class unit (grid package). Top margin.
}
  \item{margin.left}{
	Object of class unit (grid package). Left margin.
}
  \item{margin.right}{
	Object of class unit (grid package). Right margin.
}
  \item{margin.bottom}{
	Object of class unit (grid package). Bottom margin.
}
  \item{legend.display}{
	Boolean. Display or not the legend. Legend is positionned in upper right corner. 
}
  \item{legend.range}{
	Character: "full" (default) or "data". If full, color scale legend ranges from values.min to values.max. If "data", range is c(min(data),max(data)).  
}
  \item{legend.data.display}{
	Boolean. Display or not color scale legend.
}
  \item{legend.gpar}{
	Object of class gpar (grid package). Parameters to be used for color scale legend axis (cex,...).
}
  \item{legend.width}{
	Object of class unit (grid package). Width for the color scale legend.
}
  \item{legend.height}{
  	Object of class unit (grid package). Height for the color scale legend.
}
  \item{\dots}{
	Additional parameters the function may have. Not used currently.
}
}

\value{
	The function suggests a size (width, height) for the graphic returned as a vector.
	A typical usage will be to call the function a first time to get those values and call it again with an output device.
}

\author{
	Eric Lecoutre <eric.lecoutre@gmail.com>
}

\section{Colors}{
	There are several ways to specify colors used for the image zone. 
	The usual way is to have a shading from colors.groups.min to a color per group (typically the same).
	By default, a shading is indeed proposed between white (for colors.groups.min) and a same color shared by groups (red for colors.groups.max).
	The number of possible colors in the shading is determined by colors.nbreaks.
	In case one asks for distinct colors for groups, only a single value for colors.groups.min is allowed.
	By default, subgroups colors are taken from phenoData ("sampleColor" column), consequence of colors.groups being NULL.
	Colors for groups are overided by providing a vector of valid colors for this colors.groups argument.
	An additional and flexible way to determine colors is to provide a complete palette of possible colors, as a character vector of valid colors (argument colors.palette). 
	Note that in this case the argument colors.nbreaks has no effect as the number of possible values is the length of the palette.  	   
}


\examples{
	\dontrun{
library(RColorBrewer)
library(dichromat)

library(Biobase)
library(grid)
 pdf.directory=getwd()


load(file.path(getwd(),"expressionSetRma.Rda"))      #expressionSetRma


eset <- expressionSetRma[100:130,pData(phenoData(expressionSetRma))[,"sample"]%in%c(1:10,41:50)] # ARG
	##### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	exprs(eset)[1,5] <- 13.8
	exprs(eset)[10,7] <- 0.5
	eset2 <- expressionSetRma[200:250,] # ARG
	eset3 <- expressionSetRma[1000:1009,pData(phenoData(expressionSetRma))[,"sample"]%in%c(1:3,41:46)] # ARG
	eset4 <- expressionSetRma[100:230,pData(phenoData(expressionSetRma))[,"sample"]%in%c(1:20,31:50)] # ARG
	
	eset5 <- expressionSetRma[1:400,] # ARG
	
	# eset <- eset2
	
	pdf(file.path(pdf.directory,"eset.pdf"))
		size <- heatmap.expressionSet(eset,subtitle.main=" ")
	dev.off()
	pdf(file.path(pdf.directory,"eset.pdf"),width=size[1],height=size[2])
		heatmap.expressionSet(eset,subtitle.main=" ")
	dev.off()


pdf(file.path(pdf.directory,"eset2.pdf"))
	size <- heatmap.expressionSet(
			eset2,
			colors.nbreaks = 20,
			colors.pergroup=TRUE,
			legend.range="data",
			row.col.groups.display=FALSE,
			cell.gpar=gpar(lwd=0.5),
			legend.height=unit(50,"points"),
			title.just=c("center","center"),
			title.maxlines=2,
			col.groups.sep.width=unit(0,"points"),
			row.labels=featureNames(eset),
			subtitle.main="This is subtitle",
			row.order="hclust",row.groups.hclust=FALSE,
			title.gpar=gpar(cex=2),
			subtitle.gpar=gpar(cex=1.5)
	)
dev.off()
pdf(file.path(pdf.directory,"eset2.pdf"),width=size[1],height=size[2])
	size <- heatmap.expressionSet(
		eset2,
		colors.nbreaks = 20,
		colors.pergroup=TRUE,
		legend.range="data",
		row.col.groups.display=FALSE,
		cell.gpar=gpar(lwd=0.5),
		legend.height=unit(50,"points"),
		title.just=c("center","center"),
		title.maxlines=2,
		col.groups.sep.width=unit(0,"points"),
		row.labels=featureNames(eset),
		subtitle.main="This is subtitle",
		row.order="hclust",row.groups.hclust=FALSE,
		title.gpar=gpar(cex=2),
		subtitle.gpar=gpar(cex=1.5)

)
dev.off()




	
pdf(file.path(pdf.directory,"eset3.pdf"))
size <- heatmap.expressionSet(
		eset3,
		row.labels.gpar=gpar(cex=0.4,col=c(rep("red",2),rep("black",49))	), # col will correctly be a vector only if no group...
		col.labels.gpar=gpar(cex=0.6),
		colors.nbreaks = 20,
		colors.pergroup=TRUE,
		legend.range="data",
		row.col.groups.display=FALSE,
		cell.gpar=gpar(lwd=0.5),
		legend.height=unit(50,"points"),
		title.just=c("center","center"),
		title.maxlines=2,
		col.groups.sep.width=unit(0,"points"),
		row.labels=featureNames(eset),
		subtitle.main="Essai subtitle",
		row.order="hclust",row.groups.hclust=FALSE,
		interactive=FALSE
)
dev.off()

pdf(file.path(pdf.directory,"eset3.pdf"),width=size[1],height=size[2])
size <- heatmap.expressionSet(
		eset3,
		row.labels.gpar=gpar(cex=0.4,col=c(rep("red",2),rep("black",49))	), # col will correctly be a vector only if no group...
		col.labels.gpar=gpar(cex=0.6),
		colors.nbreaks = 20,
		colors.pergroup=TRUE,
		legend.range="data",
		row.col.groups.display=FALSE,
		cell.gpar=gpar(lwd=0.5),
		legend.height=unit(50,"points"),
		title.just=c("center","center"),
		title.maxlines=2,
		col.groups.sep.width=unit(0,"points"),
		row.labels=featureNames(eset),
		subtitle.main="Essai subtitle",
		row.order="hclust",row.groups.hclust=FALSE,
		interactive=FALSE
)
dev.off()



pdf(file.path(pdf.directory,"eset4.pdf"))
size <- heatmap.expressionSet(
		eset4,
		legend.range="data",
		colors.palette = dichromat(rich.colors(190)[1:128]),
		row.col.groups.display=TRUE,
		title.just=c("left","top"),
		title.maxlines=2,
		row.labels=featureNames(eset),
		subtitle.main="",
		row.order="hclust",row.groups.hclust=FALSE,
)
dev.off()

pdf(file.path(pdf.directory,"eset4.pdf"),width=size[1],height=size[2])
size <- heatmap.expressionSet(
		eset4,
		legend.range="data",
		colors.palette = dichromat(rich.colors(190)[1:128]),
		row.col.groups.display=TRUE,
		title.just=c("left","top"),
		title.maxlines=2,
		row.labels=featureNames(eset),
		subtitle.main="",
		row.order="hclust",row.groups.hclust=FALSE,
)
dev.off()

pdf(file.path(pdf.directory,"eset5.pdf"))
size <- heatmap.expressionSet(eset5,row.order="hclust",row.groups.hclust=FALSE)
dev.off()

pdf(file.path(pdf.directory,"eset5.pdf"),width=size[1],height=size[2])
heatmap.expressionSet(eset5,row.order="hclust",row.groups.hclust=FALSE)
dev.off()

	}
}
\keyword{ plot }

