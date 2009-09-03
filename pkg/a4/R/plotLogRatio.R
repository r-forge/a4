
plotLogRatio <- function(
		e,
		reference,
		within=NULL,
		across=NULL,
		nReplicatesVar=3, ## see parameters of compute function
		filename = "Rplots", # extension svg or pdf will be added according to device -- please provide full path if you want to automatically open the file
		device="svg", # svg, pdf, X11 for SVG, requires X11 and so interactive session
		
		# if col=NULL then color according to quantiles of ##TODOchange to columnsColors#test to see if quantiles over whole objects are already computed and available in expressionSetObject
		orderBy=list(rows='hclust',cols=NULL), # list with 2 parameters or NULL for nothing done at all
		# for rows: alpha, effect, hclust
		
		colorsColumns=NULL, # by default will be red  
		
		colorsColumnsBy=NULL,  # use colorsColumnsByPalette
		colorsColumnsByPalette=c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666"),# from brewer.pal(8,'accent') -- package RColorBrewer, 
		
		colorsUseMeanQuantiles = FALSE, #### If TRUE colors of bars will depends on following argument
		# test wether column quantilesColors is available in expressionSet object
		colorsMeanQuantilesPalette = c('orange','red','darkred'), # we will use the vector quantilesColors that have to be stored in expressionSet object (first use an helper function to compute quantiles on whole object)
		
		#<todo>
		colorsBarsMatrix = NULL, ## could be provided by user  -- overides colors colorsColumns and colorsColumnsBy
		colorsGenesNames = c('black'),  # colors to be used to give a color to genes names -- are recycled
		# warning: colors must be provided for genes in same order than in raw data as we don't know how they will be sorted...
		
		main = paste("log2 ratio's"),
		shortvarnames=NULL, # specify a variable from pData that will be used to display names 
		longvarnames=NULL, # if NULL, then .withingroups will be used
		gene.length = 50, 
		gene.fontsize = 6, 
		main.fontsize = 9,
		columnhead.fontsize = 8,
		mx = 1.5,
		exp.width = 1.8,
		exp.height = 0.2, 
		log2l.show = TRUE,
		log4l.show = FALSE,
		quantiles.show = FALSE,
		quantiles.compute = c(0.9),
		error.show = TRUE,
		view.psid = FALSE, 
		errorLabel = "Error bars show the pooled standard deviation", 
		closeX11=FALSE,
		openFile=FALSE,
		tooltipvalues=FALSE,
		probe2gene = TRUE,
		...){
	
	
	stopifnot(is(e,"ExpressionSet") | is(e,"ExpressionSetWithComputation"))
	
	if(!is.null(orderBy$cols)){
		if(!orderBy$cols %in% colnames(pData(e))){
			stop( "The column name by which the data has to be sorted cannot be found in the phenodata!\n")
		}
	}
	
	if (!inherits(e,"ExpressionSetWithComputation")){
		e <- try(computeLogRatio(e,reference=reference,within=within,across=across,nReplicatesVar=nReplicatesVar,...))
		if (inherits(e,'try-error')) stop('Error when doing required computation.')
		
	}  
	device <- tolower(device)
	if (!device %in% c('svg','pdf','x11','png','cairopng', 'javagd')) stop('Unknown device - use one of: svg, pdf, x11, png, javagd')
	if (device=='svg'){
		stopifnot(require(gridSVG))
		stopifnot(interactive())
	}
	if (device=='cairopng'){stopifnot(require(Cairo))}
	
	# we now work with the complete ExpressionSet object enriched with statistics (name: e)
	# there are two parts: data with averages, variances and so on and pData to handle metadata
	# we try to work as long as possible with the entire object with ExpressionSet class
	# as it automatically subset correctly rows and cols
	
	
#  browser()  
	igenes <- featureNames(e)# function in Biobase
	if (probe2gene == TRUE){
		if(is.null(featureData(e)$'Entrez ID') | is.null(featureData(e)$'Gene Symbol') | is.null(featureData(e)$'Description')){
			e <- addGeneInfo(e)
		}
		igenes.ll <- featureData(e)$'Entrez ID'
		igenes.symbol <- featureData(e)$'Gene Symbol'
		igenes.name <- featureData(e)$'Description'
		igenes.name <- paste(igenes.symbol, igenes.name, sep = " - ")
	}else{
		igenes.ll <- rep(NA,length(igenes))
		igenes.symbol <- rep(NA,length(igenes))
		igenes.name <- igenes
	}
	# if annotation is used, ensure we select only found genes
	e <- e[igenes, ]
	
	# to be able to sort short/long varnames according to columns reordering (if those arguments are provided)
	if (!is.null(shortvarnames)) names(shortvarnames)<-rownames(pData(e)[e$statistic=='diffref',])
	if (!is.null(longvarnames)) names(longvarnames)<-rownames(pData(e)[e$statistic=='diffref',])
	
	# reorder columns/rows 
	if (!is.null(orderBy)){
		if (is.list(orderBy)){
			if (!is.null(orderBy$cols)) 
				colorder <- do.call('order',
						as.list(as.data.frame(pData(e)[,orderBy$cols])))    else colorder <- 1:ncol(exprs(e))
			if (!is.null(orderBy$rows)){
				if (is.numeric(orderBy$rows))  {
					roworder <- orderBy$rows
				} else {
					if (orderBy$rows=='alpha'){  
						roworder <- order(igenes.name)
					} else if (orderBy$rows == 'effect'){
						if(sum(e$statistic == 'diffref') == 1){roworder <- order(exprs(e)[, e$statistic == 'diffref'], decreasing = TRUE)
						} else {roworder <- order(apply(exprs(e)[, e$statistic == 'diffref'], 1, FUN = function(vec){
												sum(vec)
											}), decreasing = TRUE)}
					} else if (orderBy$rows=='hclust'){
						hc=hclust(dist(exprs(e[,e$statistic=='diffref'])))
						roworder=hc$order 
					}        #D print('OK')
				}} else roworder <- 1:nrow(exprs(e))
		} else roworder <- 1:nrow(exprs(e))
		
	} else{ 
		roworder <- 1:nrow(exprs(e))
		colorder <- 1:ncol(exprs(e))
		
	}
# reorder ExpressionSet object
	e <- e[roworder,colorder]
# also order information on genes: names, colors,...
	igenes.name <- igenes.name[roworder]
	if(exists("igenes.ll")){
		igenes.ll <- igenes.ll[roworder]}
	if(exists("igenes.symbol")){
		igenes.symbol <- igenes.symbol[roworder]}
	
	
# genes names colors (recycle colors) -- sort them on roworders ! 
	colorsGenesNames <- rep(colorsGenesNames,nrow(exprs(e)))[1:nrow(exprs(e))][roworder]
	
	e.diff <- e[,e$statistic=='diffref']
	
# to be used in graph
	e.quantiles <- matrix(apply(exprs(e.diff), 2,
					function(tr) quantile(tr,probs=quantiles.compute)),nrow=length(quantiles.compute),byrow=FALSE)
	
	nc <- ncol(e.diff) # number of columns
	nr <- nrow(e.diff)
	mx <- max(abs(exprs(e.diff))) / mx
	
	
	if (is.null(shortvarnames)) shortvarnames <- pData(e.diff)$.withingroups # use adequate order as ExpressionSet object is sorted
	else shortvarnames <- shortvarnames[rownames(pData(e.diff))] #reorder 
	if (is.null( longvarnames))  longvarnames <- shortvarnames 
	else longvarnames <- longvarnames[ rownames(pData(e.diff))]      #reorder
	
	
	
	### management of colors 
	
	
# columns headers
	## order colorsColumns
	if(!is.null(colorsColumns) & all(names(colorsColumns) %in% shortvarnames)){
		colorsColumns <- colorsColumns[shortvarnames]
	}
	if (is.null(colorsColumns)|(!is.null(colorsColumnsBy))){
		if (!is.null(colorsColumnsBy)){
			uniquecolors <- as.numeric(as.factor(
							do.call('paste',c(pData(e.diff)[,colorsColumnsBy,drop=FALSE],sep='.')))) 
			if (is.null(colorsColumnsByPalette)) colorsColumns <- rainbow(max(uniquecolors))[uniquecolors]
			else colorsColumns <- rep(colorsColumnsByPalette,max(uniquecolors))[uniquecolors]
		}
		else colorsColumns <- rep('red',nc)
	}  else colorsColumns <- rep(colorsColumns,nc)[1:nc]  
	
# colorsColumns will be used for headers (by default, all red)
	
	
# bars colors 
	if (colorsUseMeanQuantiles | !is.null(colorsBarsMatrix)){
		# use colorsBarsMatrix if not NULL, else create it using the vector of colors previously built
		if (!is.null(colorsBarsMatrix)){
			# check for enough information in matrix
			stopifnot((nrow(colorsBarsMatrix)>=nr)|(ncol(colorsBarsMatrix)>=nc))
			
		}
		else {
			# use colorsMeanQuantiles recycled for each column
			# have to get genes colors from featuresData slot of ExpressionSet object
			stopifnot( "colorsQuantilesVector" %in% colnames(fData(e)))
			colorsBarsGenes <- fData(e)[,'colorsQuantilesVector']
			stopifnot(length(colorsMeanQuantilesPalette)>=max(colorsBarsGenes))
			colorsBarsMatrix <- matrix(rep(matrix(colorsMeanQuantilesPalette[colorsBarsGenes]),nc),ncol=nc)
			# we don't have to sort anymore colorsBarsMatrix as fData are already sorted
			roworder <- 1:nrow(exprs(e))
		}
	} else
	{
		colorsBarsMatrix <-  matrix(rep(colorsColumns,nr),ncol=nc,byrow=TRUE)
		# no argument specific on bars, duplicates columns headers for each bar    
		# recycle  colorsColumns vector
	}
#reorder according to genes order
	colorsBarsMatrix <- colorsBarsMatrix[roworder,,drop=FALSE]  
	
	if (rev(strsplit(filename,split="")[[1]])[4] != '.'){
		if (device=="cairopng") filename <- paste(filename,'png',sep='.')
		filename <- paste(filename,device,sep='.')
	}
	if (dirname(filename)=='.'){
		filename <- file.path(getwd(),filename) # just to ensure viewer will find it
	}
	
# open window of specified size using either X11 of PDF device
	if (device %in% c('svg', 'x11')) 
		x11(width = 3 + (nc * exp.width), height = nr * exp.height,...) else {
		if (device=='javagd') JavaGD(width = (3 + (nc * exp.width))*100, height = (nr * exp.height)*100) ## FLAG
		if (device=='pdf') pdf(filename, width = 3 + (nc * exp.width), height = nr * exp.height,...)
		if (device=='png') png(filename, width = round(100*(3 + (nc * exp.width))), height = round(100*(nr * exp.height)),...)
		if (device=='cairopng')   CairoPNG(filename, width = round(100*(3 + (nc * exp.width))), height = round(100*(nr * exp.height)),...)
	}
	##E width=min(xx,maxwidth)
	
# create main viewport
	vp1 <- viewport(x = 0, y = 0, w = 1, h = 1, just = c("left", "bottom"),
			layout = grid.layout(nrow = nr + 14, ncol = 1))
	pushViewport(vp1)
	
# draw top section of the graph (vp2)
	vp2 <- viewport(layout.pos.row = 1:6)
	pushViewport(vp2)
	
# draw gray box on top
	grid.rect(x = 0, y = 0, height = 1, width = 1, just = c("left", "bottom"),
			gp = gpar(fill = "lightgray"))
	
# add title
	grid.text(x = 0.5, y = 0.9, just = c("center", "top"), 
			gp = gpar(fontsize = main.fontsize), label = main)
	
# I reserve the first three columns of the graph for displaying the names
# of the genes, hence (i + 3)
# since I need one empty column on the right side I have a total of (nc + 4) columns
	
	
	
	
	
	for (i in 1:nc) {
		pData.i <- lapply(pData(e.diff)[i,],as.character)
		
		grid.text(x = (i + 3) / (nc + 4), y = 0.2, 
				label = shortvarnames[i],               ### color to be changed (first color)
				gp = gpar(col = colorsColumns[i], fontsize = columnhead.fontsize), just = c("center", "center"),
				name=paste("treatment_",i,sep=""))
		
		if (device=='svg'){
			pDatacol <- colnames(pData(e))[-grep('^\\.',colnames(pData(e)))]
			pDatacol <- pDatacol[pDatacol != 'statistic']
			tmp <- pData(e.diff)[,pDatacol]
			for (cn in pDatacol){
				tmp[,cn] <- paste(cn,tmp[,cn],sep='=')
			}
			details <- apply(tmp,1,paste,collapse=';')
			grid.garnish(paste("treatment_",i,sep=""),
					onmouseover=paste("highLightTreatment(evt,'",longvarnames[i],"', '",details[i],"');",sep=""),
					onmouseout="toolTip.setAttributeNS(null, 'display', 'none');" )
		}
	}
	popViewport()
	
# draw main section of the graph (vp3)
	vp3 <- viewport(layout.pos.row = 9:(nr + 9))
	pushViewport(vp3)
	for (i in 1:abs(nr / 2)) {
		grid.rect(x = 0, y = (1 / nr) * (i * 2),
				height = (1 / nr), width = 1,
				just = c("left", "centre"),
				gp = gpar(fill = colors()[63], col = 'white'))
	}
# define a vector with zeros to place multiple objects efficiently later on
	j <- rep (0, nr)
	
# main loop for placing all elements in each column
	for (i in 1:nc) {
		###E here there is a way to use a mxi rescale parameter individully for each
		### column instead of mx -- mxi to be computed using max(col)
		# make a gray box which marks 1 log2, so that scientists later can easily
		# see, how much (how many log2s) a gene was up- or down-regulated
		grid.rect(x = (i + 3) / (nc + 4), y = 0 + (1 / (nr * 2)),
				height = 1, width = (1 / (nc + 4)) / mx,
				just = c("centre", "bottom"),
				gp = gpar(fill = "lightgray", col = NULL))
		
		# draw the zero line - boxes right of this line indicate an up regulation
		# of a gene while boxes to the left indicate a down-regulation
		grid.lines(x = c((i + 3) / (nc + 4), (i + 3) / (nc + 4)),
				y = c(0 + (1 / (nr * 2)), 1 + (1 / (nr * 2))),
				gp = gpar(col = colorsColumns[i], lwd = 1.5))
		
		if ((log2l.show == TRUE) & (max(abs(exprs(e.diff)[,i]))>2)) {
			# draw the first 2 log2 line 
			grid.lines(x = c((i + 3) / (nc + 4) - (2 / (nc + 4)) / (mx * 2),
							(i + 3) / (nc + 4) - (2 / (nc + 4)) / (mx * 2)),
					y = c(0 + (1 / (nr * 2)), 1 + (1 / (nr * 2))),
					gp = gpar(col = "lightgrey", lwd = 1))
			# draw the second 2 log2 line 
			grid.lines(x = c((i + 3) / (nc + 4) + (2 / (nc + 4)) / (mx * 2),
							(i + 3) / (nc + 4) + (2 / (nc + 4)) / (mx * 2)),
					y = c(0 + (1 / (nr * 2)), 1 + (1 / (nr * 2))),
					gp = gpar(col = "lightgrey", lwd = 1))  
		}
		
		if ((log4l.show == TRUE)& (max(exprs(e.diff)[,i])>4)) {
			# draw the first 4 log2 line 
			grid.lines(x = c((i + 3) / (nc + 4) - (4 / (nc + 4)) / (mx * 2),
							(i + 3) / (nc + 4) - (4 / (nc + 4)) / (mx * 2)),
					y = c(0 + (1 / (nr * 2)), 1 + (1 / (nr * 2))),
					gp = gpar(col = "lightgrey", lwd = 1))
			# draw the second 4 log2 line 
			grid.lines(x = c((i + 3) / (nc + 4) + (4 / (nc + 4)) / (mx * 2),
							(i + 3) / (nc + 4) + (4 / (nc + 4)) / (mx * 2)),
					y = c(0 + (1 / (nr * 2)), 1 + (1 / (nr * 2))),
					gp = gpar(col = "lightgrey", lwd = 1))  
		}
		
		
		
		if (quantiles.show) {
			# draw lines for computed quantiles 
			tmp <- sapply(e.quantiles[,i],function(qq) {
						grid.lines(x = c((i + 3) / (nc + 4) + (qq / (nc + 4)) / (mx * 2),
										(i + 3) / (nc + 4) + (qq / (nc + 4)) / (mx * 2)),
								y = c(0 + (1 / (nr * 2)), 1 + (1 / (nr * 2))),
								gp = gpar(col = "lightgrey", lwd = 1))  
						
					})
			
		}
		
		# draw box according to the size of the ratio in one go for all rows
		# here comes the vector j in...
		# and it took me long to figure this one out properly...
		#E I change it to a loop to be able to add tooltips
		for (jj in 1:nrow(e.diff)){
			
#D cat('\n jj:',jj,' i:',i)
#D cat(' color:', colorsBarsMatrix[jj,i])
#D cat(' name in matrix',rownames(e.diff[jj]), ' full name:',igenes.name[jj] )
			
			grid.rect(x = (i + 3) / (nc + 4) +  
							( ( exprs(e.diff)[jj,i]  ) * ((1 / (nc + 4)) / (mx * 2)) ) / 2,
					#       y = (1 / nr)*jj, 
					y = (nr-jj+1)/nr, 
					height = (1 / (nr * 2)),
					width = abs((  exprs(e.diff)[jj,i] ) * ((1 / (nc + 4)) / (mx * 2))),
					just = c("centre", "centre"), ### change to use matrix
					
					gp = gpar(fill = colorsBarsMatrix[jj,i], col = NULL),name=paste('val',jj,i,sep='_') )
			
			if (device=='svg' & tooltipvalues){
				
				grid.garnish(paste('val',jj,i,sep='_'),
						onmouseover=paste("highLightTreatment(evt,'",round(exprs(e.diff)[jj,i],1),"', '",details[i],"');",sep=""),
						onmouseout="toolTip.setAttributeNS(null, 'display', 'none');" )
			}
			
		}  
		# draw error bar
		# only when available...
		if (error.show) {
			# retrieve column for pooled sd IF available
			#   withinData <- split(pData(e),pData(e)$.withingroups)[[i]]
			#   withinData <- withinData[withinData$statistic=='signedpooledSD',]
			errorname <- paste(
					substr(rownames(pData(e.diff))[i],1,nchar(rownames(pData(e.diff)))[i]-8)
					,'spooledSD',sep='.')
			
			
			#  if (nrow(withinData)==1){
			if (errorname %in% rownames(pData(e))){
#          psd.i <- exprs(e)[,rownames(withinData)]
				psd.i <- exprs(e)[,errorname]
				
				# horizontal line of the error bar
				# data.sd contains the data for the adjusted standard error
				# from the Dunnett's test
				grid.segments(x0 = (j + i + 3) / (nc + 4) + (exprs(e.diff)[,i]) 
								* ((1 / (nc + 4)) / (mx * 2)),
						x1 = (j + i + 3) / (nc + 4) + (exprs(e.diff)[,i]+ psd.i) 
								* ((1 / (nc + 4)) / (mx * 2)),
						#  y0 = (1 / nr) * seq(nr),
						#  y1 = (1 / nr) * seq(nr),
						
						y0 = (1 / nr) * rev(seq(nr)),
						y1 = (1 / nr) * rev(seq(nr)),
						
						
						gp = gpar(lwd = 0.5, col = "black") )
				# vertical line of the error bar
				grid.segments(x0 = (j + i + 3) / (nc + 4) +
								(exprs(e.diff)[,i] + psd.i) * ((1 / (nc + 4)) / (mx * 2)),
						x1 = (j + i + 3) / (nc + 4) + (exprs(e.diff)[,i] +psd.i) * ((1 / (nc + 4)) / (mx * 2)),
#            y0 = (1 / nr) * seq(nr) - (1 / (nr * 4)),
#            y1 = (1 / nr) * seq(nr) + (1 / (nr * 4)),
						y0 = (1 / nr) * rev(seq(nr)) - (1 / (nr * 4)),
						y1 = (1 / nr) * rev(seq(nr)) + (1 / (nr * 4)),
						
						gp = gpar(lwd = 0.5, col = "black") )
			}    
		}
	}
	
# draw the names of the genes and hyperlink them to the current annotation
# igenes.name contains a vector with the names of all the "interesting" genes
# since they are often too long, I cut them off
# igenes.ll contains the number for the entrez gene annotation
# both data are all provided by the bioconductor project
	g <- substring(igenes.name, 1, gene.length)
#  print(g)
	for (i in 1:nr) {
		grid.text(x = 0.005, 
				#y = (i/ nr), 
				y=(nr-i+1)/nr, 
				label = as.character(g[i]), name = paste("gene", i, sep = ""),
				just = c("left", "centre"), gp = gpar(fontsize = gene.fontsize, col = colorsGenesNames[i]))
		
		if (!probe2gene & device=='svg'){
			grid.hyperlink(gPath(paste("gene", i, sep = "")),  paste("http://www.ncbi.nih.gov/entrez/query.fcgi?db=gene&amp;cmd=Retrieve&amp;dopt=summary&amp;list_uids=",as.integer(igenes.ll[i]), sep = "") )
		}
	}
	popViewport()
	
	
# draw bottom section of the graph (vp4)
	vp4 <- viewport(layout.pos.row = (nr + 11):(nr + 14))
	pushViewport(vp4)
	
# draw gray box at the bottom
	grid.rect(x = 0, y = 0, height = 1, width = 1, just = c("left", "bottom"),
			gp = gpar(fill = "lightgrey"))
	
# the following commands simply draw the legend
	###<CHANGE>
	if (error.show) grid.text(x = 0.05, y = 0.5, 
				just = c("left", "bottom"), gp = gpar(fontsize = 9),label = errorLabel)
	if (colorsUseMeanQuantiles){
		grid.text(x = 0.8, y = 0.5, 
				just = c("right", "bottom"), gp = gpar(fontsize = 9),label = "Increasing quantiles ")
		k=length(colorsMeanQuantilesPalette)
		for (ik in 1:k){
			grid.rect(x = 0.8+(ik-1)*((0.9-0.8)/k),
					y = 0.5, 
					height = 0.2,
					width = (0.9-0.8)/k,
					just = c("left", "bottom"),gp = gpar(fill =colorsMeanQuantilesPalette[ik] , col = 'black'))
			
		}
	}
	
	grid.text(x = 0.05, y = 0.05, just = c("left", "bottom"),
			gp = gpar(fontsize = 4, col = "black"),
			label = paste(date(),";",R.version.string,"; Biobase version ", 
					package.version("Biobase")))
	
	popViewport()
	
# finish up the graph
# just a box around the complete graph
	grid.rect(x = 0.5, y = 0.5, width = 1, height = 1, gp = gpar(fill = NA, col = "black"))
	popViewport()
	
	if (device == 'svg') {        #E includes script for tooltips
		grid.script(file=file.path(.path.package('a4'),'etc','tooltip.script'))  
		gridToSVG(name = filename)
		if (openFile) browseURL(url=filename)
	}
	if ((device %in% c('pdf','png','cairopng')) | (device %in% c('svg','x11') & closeX11)) dev.off()
	if (device=='pdf' & openFile)  openPDF(filename)
	
	
	invisible(e)       
	
}