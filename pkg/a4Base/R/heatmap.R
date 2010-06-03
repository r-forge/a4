#<TODO> proposed adapted cex for text display according to size (ncol/nrow) of exprs object 
# --> better adapt proposed height/width! 


# FUNCTIONS
splitTitle <- function(vec,cp=30,maxlines=4){
	a1= sapply(vec, strsplit," ")
	a2= lapply(a1,FUN=function(el) cumsum(nchar(el)))
	a3= lapply(a2,FUN=function(vec) floor(vec/cp)+1)
	nlignes=sapply(a3,max)
	newlab <- vector(mode="list",length=length(a3))
	for (i in 1:length(newlab)){
	    if (nlignes[i]>1){
	      a4 <-c(lapply(split(a1[i][[1]],as.factor(a3[i][[1]])),paste,collapse=" "))
	      if (length(a4)>maxlines) {
	        if (nchar(a4[maxlines])>cp) {
					a4[maxlines] <- paste(substr(a4[maxlines],1,cp-3),"/...",sep="")
				} else a4[maxlines] <- paste(a4[maxlines],"/...",sep="")
	        a4 <- a4[1:maxlines]
	       }
	      a5 <- do.call("paste",c(a4,sep="\n"))
	      newlab[i] <- a5
	    }
	    else newlab[i] <- vec[i]
	}
	
	tmplab <- unlist(newlab)
	return(tmplab)
}


# slighlty adapted from P. Murrell "R Graphics", pages 233
makeImageRect <- function(nrow, ncol,cols, byrow,gp=list(col=NULL),
			force.y=NULL,force.height=NULL, just=c("right", "top")) {
  xx <- (1:ncol)/ncol   
  yy <- (1:nrow)/nrow
  
  if (byrow) {
    right <- rep(xx, nrow)
    top <- rep(yy, each=ncol)
  } else {
    right <- rep(xx, each=nrow)
    top <- rep(yy, ncol)
  }  
  width=1/ncol
  height=1/nrow
  if (!is.null(force.y)) top <- force.y
  if (!is.null(force.height)) height <- force.height
  
  rectGrob(x=right, y=top, 
           width=width, height=height, 
           just=just, 
           gp=gpar(fill=cols,col=NULL),
           name="image")
}

imageGrob <- function(nrow, ncol, cols, byrow=TRUE,
                       name=NULL, gp=NULL, vp=NULL,...) { 
  igt <- gTree(nrow=nrow, ncol=ncol, 
               cols=cols, byrow=byrow,
               children=gList(makeImageRect(nrow, ncol, 
                                            cols, byrow,...)),
               gp=gp, name=name, vp=vp, 
               cl="imageGrob") 
  igt
}

grid.imageGrob <- function(...) {
  igt <- imageGrob(...)
  grid.draw(igt)
}


heatmap.expressionSet <- function(
	eset,
	col.groups = pData(phenoData(eset))[,"subGroup"],
	col.orderBy = order(pData(phenoData(eset))[,"subGroup"]), 
	col.groups.sep.width=unit(8,"points"),
	col.labels = sampleNames(eset),
	col.labels.sep.width=unit(10,"points"),
	col.labels.gpar = gpar(cex=1), 
	col.labels.max.nchar = 20,
	colors.pergroup=FALSE,
	colors.groups = NULL,
	colors.groups.min  = rgb(1,1,1),
	colors.max = rgb(1,0,0),
	colors.min = rgb(1,1,1),
	colors.nbreaks = 128,
	colors.palette = NULL,
	cell.gpar = gpar(lty=0),
	row.groups.sep.height=unit(15,"points"),
	row.labels.sep.height=unit(10,"points"),
	row.col.groups.display = ifelse(length(unique(col.groups))>1,TRUE, FALSE),
	row.col.groups.display.height = unit(6,"points"),
	row.labels.gpar = gpar(cex=1,col="black"),
	row.labels.max.nchar = 45,
	row.labels = list("SYMBOL","GENENAME"),#row.labels <- featureNames(eset)
	row.labels.sep=" - ",
	row.groups = rep(1,nrow(exprs(eset))), 
	row.order="none", # hclust, integer bvector
	row.groups.hclust = FALSE, # if row.order="hclust",possibility to split groups
	row.groups.hclust.n = 4,
	distfun       = dist,
	hclustfun     = function(d) { hclust(d, method = "ward") },
	values.min=0,
	values.max=16,
	title.gpar = gpar(cex=1.4),
	title.main="This is the title  possibly being very long - it will be splited on several lines or even displayed with dots at the end -- see there (does it work? addendum)",
	title.just=c("right","top"), # ! has to be a vector of size 2
	title.maxlines=4,
	title.cutpoint=40,
	subtitle.gpar=gpar(cex=1),
	subtitle.main="This is subtitle",
	subtitle.maxlines=4,
	subtitle.just=title.just, # ! has to be a vector of size 2,
	subtitle.cutpoint=40,
	margin.top = unit(2,"lines"),
	margin.left = unit(2,"lines"),
	margin.right = unit(2,"lines"),
	margin.bottom = unit(2,"lines"),
	legend.display=TRUE,
	legend.range="full", # or: data
	legend.data.display=ifelse(legend.range=="full",TRUE,FALSE),
	legend.gpar = gpar(cex=1),
	legend.width = unit(250,"points"),
	legend.height  = unit(40,"points")	
	,...
	){

		
		stopifnot(require(Biobase))
		stopifnot(require(grid))
		
		### REORDER EXPRESSION SET
		################################################################################
		#row.order=match.arg(row.order,c("hclust","none"))
		
		if (is.numeric(row.order))rowInd <- row.order else{
			if (row.order=="none"){
			    rowInd=1:nrow(exprs(eset))
			} else if (row.order=="hclust"){
			    #from heatmap2 -- to order rows according to foo
				x=exprs(eset)
			    Rowv <- rowMeans(x, na.rm = TRUE)
			    hcr <- hclustfun(distfun(Rowv))
			    ddr <- as.dendrogram(hcr)
				if (row.groups.hclust){
					row.groups <- cutree(hcr,row.groups.hclust.n) #cutree(hcr,3) -- for groups
					#ddr <- reorder(ddr, Rowv)
					rowInd <- order(row.groups,Rowv)#,order.dendrogram(ddr))
					
				}else {
					ddr <- reorder(ddr, Rowv)
					rowInd <- order.dendrogram(ddr)
				}
			} else stop("Row order asked not implemented -- use one of: numeric vector, 'none','hclust'") 
		}
		
		eset <- eset[rowInd,col.orderBy]
		
		################################################################################
		
		# tmp 
		data <- exprs(eset)
		col.groups.n=length(unique(col.groups))
		row.groups.n=length(unique(row.groups))
		col.groups.table <- table(col.groups)
		row.groups.table <- table(row.groups)
		
		                                    
		
		### Preparation of colors
		################################################################################
		
		if (is.null(colors.groups))
		{
			#print("!colors.groups")
			tmp=aggregate(pData(phenoData(eset))[,"sampleColor"],by=list(g=pData(phenoData(eset))[,"subGroup"]),FUN=unique)
			colors.groups <- as.character(tmp[order(tmp$g),"x"])	
		}
		if (is.null(colors.palette)){
			#print("!colors.palette")
			# to ensure we can always use colors
			colors.max <-rep(colors.max,col.groups.n)[1:col.groups.n]
			colors.min <-rep(colors.min,col.groups.n)[1:col.groups.n]
			  if (colors.pergroup){
				  #print("There1")
			      colors.2break <- cbind(
			          t(col2rgb(colors.groups,alpha=TRUE)),
			          t(col2rgb(rep(colors.groups.min,col.groups.n),alpha=TRUE)))
			  colnames(colors.2break) <- NULL
			        
					  colors.2break <- as.data.frame(t(colors.2break))
			        colors.groups.shade=lapply(colors.2break,FUN=function(vec.minmax){
			          colorpanel(n=colors.nbreaks,
			            high= do.call("rgb",c(as.list(vec.minmax[1:4]),maxColorValue=255)),
			            low=do.call("rgb",c(as.list(vec.minmax[5:8]),maxColorValue=255)))
			            })
			   } else
			  {
					#print("There2")
			      colors.2break <- cbind(
			          t(col2rgb(colors.max,alpha=TRUE)),
			          t(col2rgb(colors.min,alpha=TRUE)))
			          colnames(colors.2break) <- NULL
			          colors.2break <- as.data.frame(t(colors.2break))
			        colors.groups.shade=lapply(colors.2break,FUN=function(vec.minmax){
			          colorpanel(n=colors.nbreaks,
			            high= do.call("rgb",c(as.list(vec.minmax[1:4]),maxColorValue=255)),
			            low=do.call("rgb",c(as.list(vec.minmax[5:8]),maxColorValue=255)))
			            })
					}
			  } else{
				  
				  tmp <- list()
				 for (i in 1:col.groups.n){
						tmp[[i]] <- colors.palette
					  
				  }
				  #print("palette!")
				  colors.groups.shade<-tmp
				  colors.nbreaks <- length(colors.palette)
			  }	  

		# break data into pieces and assign colors
		breaks = seq(from=values.min,to=values.max,length.out=colors.nbreaks)
		
		min.breaks <- min(breaks)
		max.breaks <- max(breaks)
		data[data < min.breaks] <- min.breaks
		data[data > max.breaks] <- max.breaks
		# from image.default --> to cut into classes
		data.breaks <- matrix(.C("bincode", as.double(data), length(data), as.double(breaks),
		            length(breaks), code = integer(length(data)), (TRUE),
		            (TRUE), nok = TRUE, NAOK = TRUE, DUP = FALSE, PACKAGE = "base")$code -
		            1 ,ncol=ncol(data), nrow=nrow(data),byrow=FALSE)+1
		
		
		#print(breaks)			
		data.breaks.groups0 <- split(as.data.frame(t(data.breaks)),as.factor(col.groups))
		data.breaks.groups1 <- list()
		for (i in 1:length(data.breaks.groups0)){
			data.breaks.groups1[[i]] <- as.data.frame(t(data.breaks.groups0[[i]]))
			
		}
		#data.breaks.groups1 <- lapply(),FUN=function(group) as.data.frame(t(group)))
		data.breaks.groups <- list()
		for (j in 1:length(data.breaks.groups1)){
			data.breaks.groups[[j]] <- split(data.breaks.groups1[[j]],as.factor(row.groups))
			
		}
	#data.breaks.groups <- lapply(data.breaks.groups1,FUN=function(group) split(group,as.factor(row.groups)))
		
		colors.breaks.groups <- data.breaks.groups
		for (colGroup in 1:length(data.breaks.groups)) {
		  for (rowGroup in names(data.breaks.groups[[colGroup]])){ 
		      colors.breaks.groups[[colGroup]][[rowGroup]] <- 
				matrix(colors.groups.shade[[colGroup]]
					[as.matrix(
						data.breaks.groups[[colGroup]][[rowGroup]])],
						ncol=ncol(data.breaks.groups[[colGroup]][[rowGroup]]))
		  	}
		}
		
		##############################################################################
		### LABELS   -- from heatMapIntensities
		### Note: xlab and ylab changed to row.labels and col.labels
		##############################################################################
		
		### clean title
		
		title.toplot <-splitTitle(title.main,cp=title.cutpoint,maxlines=title.maxlines)
		title.toplot.grob <- textGrob(title.toplot,gp=title.gpar)
		
		if (subtitle.main=="") {subtitle.main=" "}
		subtitle.toplot <- splitTitle(subtitle.main,cp=subtitle.cutpoint,maxlines=subtitle.maxlines)
		subtitle.toplot.grob <- textGrob(subtitle.toplot,gp=subtitle.gpar)
		
		  # retrieve row labels from ExpressionSet object
			if (!is.list(row.labels)) {
					row.labels <- row.labels[rowInd] # reorder names
			}	else {
				tmpRowLabels <- unlist(row.labels)###
				tmpRowLabels <- tmpRowLabels[tmpRowLabels %in% colnames(pData(featureData(eset)))]
				if (length(tmpRowLabels)==0){
					stop("Columns of pData selected for row names do not exist in expressionSet")
				} else {
					index <- rownames(data)
					tmp<-as.matrix(pData(featureData(eset))[index,tmpRowLabels,drop=FALSE])
					tmp[is.na(tmp)] <- ""
					row.labels2 <- apply(tmp,1,paste,collapse=row.labels.sep)
					row.labels2[gsub(" ", "", row.labels2)==""] <- names(row.labels2[gsub(" ", "", row.labels2)==""]) 
					row.labels <- row.labels2			
				}
			}
			if (!is.list(col.labels)) {
				col.labels <- col.labels[col.orderBy] # reorder names
			}	else {
				
				tmpColLabels <- unlist(col.labels)###
				tmpColLabels <- tmpColLabels[tmpColLabels %in% colnames(pData(phenoData(eset)))]
				if (length(tmpColLabels)==0){
					stop("Columns of pData selected for samples do not exist in expressionSet")
				} else {
					index <- colnames(data)
					tmp<-as.matrix(pData(phenoData(eset))[index,tmpColLabels,drop=FALSE])
					tmp[is.na(tmp)] <- ""
					col.labels2 <- apply(tmp,1,paste,collapse=' ')
					col.labels2[gsub(" ", "", col.labels2)==""] <- names(col.labels2[gsub(" ", "", col.labels2)==""]) 
					col.labels <- col.labels2			
				}
			}
			
		   # clean labels: cut them if too long
		  row.labels.n=length(row.labels)
		  row.labels.toplot <- sapply(row.labels,FUN=function(str){
		      if (nchar(str) >row.labels.max.nchar) {
		      paste(  substr(str,1,row.labels.max.nchar),"[...]",sep="")
		      } else str
		    })
		  row.labels.str.max <- max(nchar(row.labels.toplot))
		  row.labels.str.longer <-(row.labels.toplot[which(nchar(row.labels.toplot)==row.labels.str.max)])[1]
		  row.labels.toplot.groups <- split(row.labels.toplot,as.factor(row.groups))
		
		  col.labels.n=length(col.labels)
		  col.labels.toplot <- sapply(col.labels,FUN=function(str){
		      if (nchar(str) >col.labels.max.nchar) {
		      paste(  substr(str,1,col.labels.max.nchar),"[...]",sep="")
		      } else str
		    })
		  col.labels.str.max <- max(nchar(col.labels.toplot))
		  col.labels.str.longer <-(col.labels.toplot[which(nchar(col.labels.toplot)==col.labels.str.max)])[1]
		  col.labels.toplot.groups <- split(col.labels.toplot,as.factor(col.groups))
		
		################################################################################
		
		
		#### PREPARE LAYOUT FOR BLOCS OF IMAGES
		################################################################################
		
		
		tmp.forgroups <- if(col.groups.n>1){
		    values <- c()
		    for (i in 2:col.groups.n){
		      values <- c(values,col.groups.sep.width[[1]],col.groups.table[i])
		    }
		    values  
		  } else NULL
		tmp.forgroups.unit <- if(col.groups.n>1){rep(c(attr(col.groups.sep.width,"unit"),"null"),col.groups.n-1)} else NULL
		tmp.forgroups.data.nNULL <- if(col.groups.n>1){2*(col.groups.n-1)} else 0
		
		col.units.values <- c(col.groups.table[1],tmp.forgroups,col.labels.sep.width[1],1)
		col.units.units  <- c("null",tmp.forgroups.unit,attr(col.labels.sep.width,"unit"),"grobwidth")
		col.units.data <- c( vector(mode="list",length=1+tmp.forgroups.data.nNULL+1),list(textGrob(label=row.labels.str.longer,gp=row.labels.gpar)))
		
		grid.layout.col.widths=unit(
		    col.units.values,
		    col.units.units,
		    col.units.data)
		
		#tmp.forgroups<- if(row.groups.n>1){rep(c(row.groups.sep.height[1],1),row.groups.n-1)} else NULL
		#tmp.forgroups.unit <- if(row.groups.n>1){rep(c(attr(row.groups.sep.height,"unit"),"null"),row.groups.n-1)} else NULL
		#tmp.forgroups.data.nNULL <- if(row.groups.n>1){2*(row.groups.n-1)} else 0
		tmp.forgroups <- if(row.groups.n>1){
					values <- c()
					for (i in 2:row.groups.n){
						values <- c(values,row.groups.sep.height[[1]],row.groups.table[i])
					}
					values  
				} else NULL
		tmp.forgroups.unit <- if(row.groups.n>1){rep(c(attr(row.groups.sep.height,"unit"),"null"),row.groups.n-1)} else NULL
		tmp.forgroups.data.nNULL <- if(row.groups.n>1){2*(row.groups.n-1)} else 0
		
		
		row.units.values <- c(
				row.groups.table[1],
				tmp.forgroups,
				if(row.col.groups.display){c(6,row.col.groups.display.height[1])},
				row.labels.sep.height[1],
				1)
		row.units.units <- c(
				"null",
				tmp.forgroups.unit,
				if(row.col.groups.display)c("points",attr(row.col.groups.display.height,"unit")),
				attr(row.labels.sep.height,"unit"),
				"grobheight")
		row.units.data <- c( vector(mode="list",length=1+tmp.forgroups.data.nNULL+1+{ifelse(row.col.groups.display,2,0)}),list(textGrob(col.labels.str.longer,gp=col.labels.gpar,rot=90)))
		
		grid.layout.row.heights=unit(
		    row.units.values,
		    row.units.units,
		    row.units.data  )
		
		heatmap.grid.nrow = 1+2*(row.groups.n-1)+1+2*row.col.groups.display+1
		heatmap.grid.ncol=1+2*(col.groups.n-1)+2
		
		heatmap.grid = grid.layout(  
		  nrow=heatmap.grid.nrow,
		  ncol=heatmap.grid.ncol,
		  width=grid.layout.col.widths,
		  heights=grid.layout.row.heights)
		
		
		# prepare indexes of blocs where to plot "images" of sub-groups
		# take into account the gap between viewports in layout
		
		#  blocs.col.index=1:col.groups.n
		#  if  (col.groups.n>1) {
		#    for (i in 1:(col.groups.n-1)){ blocs.col.index <- blocs.col.index+(1*blocs.col.index>i)}
		#  # todo test with more groups
		#  }
		
		  blocs.col.index<- seq(1,2*col.groups.n,by=2)
		  blocs.row.index<- seq(1,2*row.groups.n,by=2)
		  
		  blocs.all <- expand.grid(y=blocs.row.index,x=blocs.col.index)
		  tmp <- expand.grid(yind=1:row.groups.n,xind=1:col.groups.n)
		  blocs.all <- cbind(bloc=1:nrow(blocs.all),blocs.all,tmp)
		
		
		
		# grid.show.layout(heatmap.grid)
		# print(blocs.all)
		
		
		#### PREPARE TOP LAYOUT 
		################################################################################
		
		#if (interactive) {
			grid.newpage()
		#	}
		
		title.legend.height.points <- max(
			convertUnit(legend.height+unit(40,"points"),"points")[[1]], 
			convertUnit(unit(1,"grobheight",list(title.toplot.grob)),"points")[[1]]
				+ convertUnit(unit(1,"grobheight",list(subtitle.toplot.grob)),"points")[[1]]
	)
		
		top.vp.layout <- grid.layout(
		  nrow=5, ncol=3,
		 widths = unit(
		    c(margin.left[[1]], 1,margin.right[[1]]),
		    c(attr(margin.left,"unit"), "null",attr(margin.right,"unit"))),
		 heights = unit(
		    c(margin.top[[1]],  title.legend.height.points ,1,1,margin.bottom[[1]]),
			c(attr(margin.top,"unit"), "points","lines","null", attr(margin.bottom,"unit")), 
			
			))
		
		top.vp <- viewport(layout=top.vp.layout)
		# grid.show.layout(top.vp.layout)
		
		pushViewport(top.vp)
		
		
		# push viewport with the layout for images
		heatmap.plot.zone <- viewport(
		  layout=heatmap.grid,
		  layout.pos.col=2,layout.pos.row=4)
		pushViewport(heatmap.plot.zone)
		
		
		# draw blocs
		for (ibloc in blocs.all[,"bloc"]){
			ibloc.colors <- colors.breaks.groups[[blocs.all[ibloc,"xind"]]][[blocs.all[ibloc,"yind"]]]
			#grid.rect(vp=viewport(layout.pos.row=blocs.all[ibloc,"y"],layout.pos.col=blocs.all[ibloc,"x"]))
			ibloc.colors  <- ibloc.colors[rev(1:nrow(ibloc.colors)),,drop=FALSE]  
			grid.imageGrob(
					nrow(ibloc.colors),
					ncol(ibloc.colors),
					as.vector(t(ibloc.colors)),
					vp=viewport(layout.pos.row=blocs.all[ibloc,"y"],layout.pos.col=blocs.all[ibloc,"x"]),
					gp=cell.gpar,byrow=TRUE)
			# labels
			tmplabels <- col.labels.toplot.groups[[blocs.all[ibloc,"xind"]]]
			grid.text(tmplabels,
					x=(0:(length(tmplabels)-1))/length(tmplabels)+1/(2*length(tmplabels)),
					just=c("center","center"),gp=col.labels.gpar,rot=90,draw = TRUE, 
					vp = viewport(layout.pos.row=heatmap.grid.nrow,layout.pos.col=blocs.all[ibloc,"x"]))
			
			tmplabels <- row.labels.toplot.groups[[blocs.all[ibloc,"yind"]]]
			
			grid.text(tmplabels,
					x=0,
					#y=(((length(tmplabels)):1)-0.5)/length(tmplabels),#- 1/(2*length(tmplabels))
					y=seq(1,1/(2*length(tmplabels)),by=-1/(length(tmplabels)))-1/(2*length(tmplabels))
					,
					just=c("left","center"),gp=row.labels.gpar,draw = TRUE, 
					
					vp = viewport(
							
							layout.pos.row=	blocs.all[ibloc,"y"],
							layout.pos.col=heatmap.grid.ncol))
		}  
		
		if (row.col.groups.display){
			for (igroup in 1:col.groups.n){
				#  	print(igroup)
				grid.rect(gp=gpar(fill=colors.groups[igroup],lty=0),
						vp=viewport(layout.pos.row=heatmap.grid.nrow-2,layout.pos.col=blocs.col.index[igroup]))
			}
		}
			# come back to main viewport
		popViewport()
		
		## title #
		
		# push viewport with the layout for images
		title.vp.layout <-grid.layout(
				nrow=1,
				ncol=5,
				width=unit(
							c(1,1,30,legend.width[1],1),
							c("null","grobwidth","points",attr(legend.width,"unit"),"null"),
							c(vector(mode="list",length=1),list(title.toplot.grob),vector(mode="list",length=3))))
		
		# grid.show.layout(title.vp.layout)
		
		title.plot.zone <-
			viewport(
				layout=title.vp.layout,
				layout.pos.col=2,layout.pos.row=2)
		pushViewport(title.plot.zone)
		# grid.rect(gp=gpar(col="red"))

		
		
		title.vp <- viewport(layout.pos.row=1,
					layout.pos.col=c(1,2),
					layout=grid.layout(3,1,
						heights=unit.c(
							unit(1,"grobheight",list(title.toplot.grob))
							,unit(1,"lines")
							,	unit(1,"grobheight",list(subtitle.toplot.grob))
								)))
			
		pushViewport(title.vp)
			pushViewport(viewport(layout.pos.row=1,layout.pos.col=1))
				# grid.rect(gp=gpar(col="blue"))
				tmp.just.x <- match(title.just[1],c("center","left","right"))
				tmp.x <- c(0.5,0,1)
				
				tmp.just.y <- match(title.just[2],c("center","bottom","top"))
				tmp.y <- c(0.5,0,1)
				
				grid.text(title.toplot,x=tmp.x[tmp.just.x],y=tmp.y[tmp.just.y],just=title.just,gp=title.gpar)
			popViewport()
			pushViewport(viewport(layout.pos.row=3,layout.pos.col=1))
				# grid.rect(gp=gpar(col="green"))
				tmp.just.x <- match(subtitle.just[1],c("center","left","right"))
				tmp.x <- c(0.5,0,1)
				
				tmp.just.y <- match(subtitle.just[2],c("center","bottom","top"))
				tmp.y <- c(0.5,0,1)
				
				grid.text(subtitle.toplot,x=tmp.x[tmp.just.x],y=tmp.y[tmp.just.y],
							just=subtitle.just,gp=subtitle.gpar)
			popViewport()
			
			
		popViewport()
		if (legend.display){
		
				legend.vp <- viewport(layout.pos.row=1,layout.pos.col=4,
							layout=grid.layout(2,1,heights=unit.c(legend.height+unit(5,"points"),unit(1,"null"))))
				pushViewport(legend.vp)
				legend.draw.vp <- viewport(layout.pos.row=1,layout.pos.col=1)
				pushViewport(legend.draw.vp)
		#	grid.rect(gp=gpar(col="green"))
			
			if (legend.range=="full"){
				# display full range
				tmpNcol <- length(colors.groups.shade[[1]])
		
				grid.imageGrob(nrow=1,ncol=tmpNcol,cols=colors.groups.shade[[1]],
						byrow=TRUE,gp=cell.gpar,force.y=1,force.height=legend.height)
				grid.rect(x=0.5,y=1,height=legend.height,width=unit(1,"npc"),just=c("center","top"))

				if (legend.data.display){
					rr=(range(data)-values.min)/values.max
		#			grid.lines(x=c(rr[1],rr[1]),y=unit.c(unit(1,"npc")-convertUnit(legend.height,"npc"),unit(1,"npc")),gp=gpar(lwd=2))
		#			grid.lines(x=c(rr[2],rr[2]),y=unit.c(unit(1,"npc")-convertUnit(legend.height,"npc"),unit(1,"npc")),gp=gpar(lwd=2))
					
					grid.lines(x=c(rr[1],rr[1]),y=unit.c(unit(1,"native"),unit(0,"native")+unit(5,"points")),gp=gpar(lwd=2))
					grid.lines(x=c(rr[2],rr[2]),y=unit.c(unit(1,"native"),unit(0,"native")+unit(5,"points")),gp=gpar(lwd=2))
					
				}
				x=c(values.min,pretty(c(values.min,as.vector(data),values.max)))
				xx <- pretty(x)
				xx[xx<values.min]<- values.min
				xx[xx>values.max]<- values.max
				xx <- unique(c(values.min,signif(range(data)[1],3),xx,signif(range(data)[2],3),values.max))
				xx.where <- (xx-min(xx))/max(xx)
				grid.xaxis(at=xx.where,label=xx,gp=legend.gpar)
			} else {
				# display only range of data
				tmpseq=seq(min(data.breaks),max(data.breaks),by=1)
				tmpNcol=length(tmpseq)
						
				grid.imageGrob(nrow=1,ncol=tmpNcol,cols=colors.groups.shade[[1]][tmpseq],
						byrow=TRUE,gp=cell.gpar,force.y=1,force.height=legend.height)
				grid.rect(x=0.5,y=1,height=legend.height,width=unit(1,"npc"),just=c("center","top"))
				xx=pretty(breaks[tmpseq])
				xx=xx[-c(1:2,length(xx),length(xx)-1)]
				xx <- signif(c(min(data),xx,max(data)),digits=3)
				x.where <- (xx-breaks[min(data.breaks)])/(breaks[max(data.breaks)+1]-breaks[min(data.breaks)]) 
				grid.xaxis(at=x.where,label=xx,gp=legend.gpar)
			}
			}
			popViewport()		
		popViewport()
	popViewport()

	
	# compute propose width/height for plot
	
#row.groups.sep.height=unit(15,"points"),
#row.labels.sep.height=unit(10,"points"),
#row.col.groups.display = ifelse(length(unique(col.groups))>1,TRUE, FALSE),
#row.col.groups.display.height = unit(6,"points"),

	proposed.height=margin.top+margin.bottom
	if (title.main!="") {proposed.height = proposed.height+unit(1,"strheight",list(title.toplot))}
	if (subtitle.main!="") {proposed.height = proposed.height+unit(1,"strheight",list(subtitle.toplot))}
	if (row.col.groups.display)	{proposed.height = proposed.height+row.col.groups.display.height}
	proposed.height = proposed.height+unit(1,"lines") #between title and main zone

	row.proposed.units.values <- row.units.values
	row.proposed.units.values[sapply(row.units.units,FUN=function(el){el})=="null"] <- 
			row.proposed.units.values[sapply(row.units.units,FUN=function(el){el})=="null"]*30
	row.proposed.units.units <- row.units.units
	row.proposed.units.units[sapply(row.units.units,FUN=function(el){el})=="null"] <- "points"
	row.proposed.units.data <-row.units.data 
	row.proposed.units <- unit(row.proposed.units.values,row.proposed.units.units,row.proposed.units.data)
	row.proposed.heights <-row.proposed.units[1]
	for (i in 2:length(row.proposed.units)){row.proposed.heights=row.proposed.heights+row.proposed.units[i]}
	proposed.height = proposed.height+row.proposed.heights

	

	proposed.width.1=margin.left+margin.right
	proposed.width.1=proposed.width.1+max(unit(1,"grobwidth",title.toplot.grob),unit(1,"grobwidth",subtitle.toplot.grob))
	proposed.width.1=proposed.width.1+unit(30,"points")
	proposed.width.1=proposed.width.1+legend.width #scale legend

	
	proposed.width.2=margin.left+margin.right

	col.proposed.units.values <- col.units.values
	col.proposed.units.values[sapply(col.units.units,FUN=function(el){el})=="null"] <- 
			col.proposed.units.values[sapply(col.units.units,FUN=function(el){el})=="null"]*25
	col.proposed.units.units <- col.units.units
	col.proposed.units.units[sapply(col.units.units,FUN=function(el){el})=="null"] <- "points"
	col.proposed.units.data <-col.units.data 
	col.proposed.units <- unit(col.proposed.units.values,col.proposed.units.units,col.proposed.units.data)
	col.proposed.width <-col.proposed.units[1]
	for (i in 2:length(col.proposed.units)){col.proposed.width=col.proposed.width+col.proposed.units[i]}
	proposed.width.2 = proposed.width.2+col.proposed.width
	
	proposed.width=max(proposed.width.1, proposed.width.2)

	proposed.height.points <- convertUnit(proposed.height,"inches",valueOnly=TRUE)
	proposed.width.points <- convertUnit(proposed.width,"inches",valueOnly=TRUE)
	
	return(c(proposed.width.points,proposed.height.points))
}

