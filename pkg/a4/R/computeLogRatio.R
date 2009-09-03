# e: expressionSet object
# reference$var: variable that for which computations (averages, diff to control,...) will be computed and displayed displayed in graph
# reference$level: specific level of var from which differences are computed
# across: vector of variables that make combinations of grouping for computations
# var and across variables MUST be present in pData(ExpressionSetObject)
# within: reference$var and across could have same combinations in different subgroups of columns such as those that come from celllines - we first within data, make computations per within and merge back at the end ; this is a vector of withinting variables

## Note: if there  are additional pheno data than the ones used in within,var,across their content could be lost in aggregated pData used to store statistics

## <DEBUG>
# data(e)
# reference=list(var='dose',level=10)
#  within=c('cellline')
# across=NULL
# nReplicatesVar=3

## </DEBUG>


'computeLogRatio' <- function(e,reference,within=NULL,across=NULL,nReplicatesVar=3,...){

  stopifnot(is(e,'ExpressionSet'))
  stopifnot(colnames(e)==sampleNames(e))
  if (!all(c(reference$var,reference$across)%in% colnames(pData(e)))) 
  {
    stop(paste('\n reference var or any variable in across not present in phenotype data',' (',
        paste(c(var,across),collapse=' '),')\n',sep=''))
  }
  if (! reference$level %in% unique(pData(e)[,reference$var])) 
  {
    stop(paste('\n reference level is not define for reference$var'))
  }

 # if(is.null(within)) {
 #   pData(e)<-cbind(pData(e),'.tmpwithin'=1)
 #   
 # }
  pData <- pData(e)
  pData <- pData[do.call('order',as.list(pData[,c(within,reference$var,across),drop=FALSE])),]
  # here we use names with . just in case phenoData already contains such a variable
  # adding reference level of var in pData to use it in further computations (make the code more readable)

  pData[,'.reference'] <- pData[,reference$var]==reference$level
  # careful if NA is one of the possible levels for var... then our variable reference would have some NA
  if (is.na(reference$level)) pData$.reference[is.na(pData$.reference)] <- TRUE else
    pData$.reference[is.na(pData$.reference)] <- FALSE 
  
  # adding combinations of within reference$var & across for further combinations
  pData$.groups <- do.call('paste',
    c(as.list(pData[,c(reference$var,across),drop=FALSE]),sep='.'))
  


  # adding replicates with the selected parameters
  pData$.replicates <- replicates(do.call('paste',
    c(as.list(pData[,c(within,reference$var,across),drop=FALSE],sep='.'))))
  
  # preparing withinted data, (could be more than one within variables)
  if (is.null(within)){
    bywithin <- rep(1,nrow(pData))
  } else {
    bywithin <- do.call('paste', 
      c(as.list(pData[,within,drop=FALSE]),sep='.'))
      pData$.within <- bywithin
  }

  if (is.null(within)) pData$.withingroups <- paste(pData$.groups,sep='.') else pData$.withingroups <- paste(bywithin,pData$.groups,sep='.')
  pDatawithin <- split(pData,bywithin) 

  # check there is at least one reference in each within groups
  check <- sapply(pDatawithin, function(block) any(block$.reference))
  if (any(!check)) stop('Error in computation: within variables define groups for which there are not any reference to make compare with.')

  # check if defined reference level and across are compatible: level must define
  # groups as little as possible - in other words, we can't have two groups 'across' 
  # that are reference level for a certain within
   check <- any(
      sapply(split(pData,pData$.within),
        function(biggroup) length(unique(biggroup[biggroup$.reference,'.groups']))>1))

   if (check) stop('Incompatibility between reference level defined by var and across variables. \nReference level must constitute a group as little as possible.\nIt appears here one across variable splits the reference level into several subgroups.\nPlease modify your call.')

  exprs.new <- c()
  pData.new <- c()

  for (iwithin in seq(along=pDatawithin)){
    # cat('\n***DEBUG - iwithin: ',iwithin)
    #E iwithin <- 1 # iwithin <- 2
    iwithinname <- names(pDatawithin)[iwithin]
    ipData <- pDatawithin[[iwithin]]
    ipData <- ipData[order(ipData$.withingroups),] # necessary so that further tapply do have same order than within (alphabetic order)
 
 
    # prepare pData that will be used for computed aggregated statistics
    ireference <- tapply(ipData$.reference,ipData$.groups,function(ref)sum(ref)>0)   
    ireplicates <- tapply(ipData$.replicates,ipData$.groups,max)   

    statpData <- ipData[!duplicated(ipData[,c(within,reference$var,across),drop=FALSE]),,drop=FALSE]
      # remove columns that were added before as they will be replaced
    statpData <- statpData[,-which(colnames(statpData)=='.reference')]
    statpData <- statpData[,-which(colnames(statpData)=='.replicates')]
    if ('.oldcolnames' %in% colnames(statpData)){
      statpData <- statpData[,-which(colnames(statpData)=='.oldcolnames')]
    }
    statpData <- cbind(statpData,.reference=ireference,.replicates=ireplicates)

    # withinting data in groups of across (within within) to compute statistics per groups in a row
    igroups <- split(ipData,ipData$.groups) 

    # prepare root names that will be used to store statistics

    irootnames <-  do.call('paste',
      c(as.list(statpData[,c(within,reference$var,across),drop=FALSE]),sep='.'))
    
#    irootnames <-  do.call('paste',
#     c(as.list(statpData[,c(within,across),drop=FALSE]),sep='.'))
    # here begin the computations 
    # each time we also prepare adequate pData
  
    ## computing averages    
    exprs.mean <- do.call('cbind',
        lapply(igroups,function(group) 
            apply(exprs(e)[,rownames(group)],1,mean)))
    colnames(exprs.mean) <- paste(irootnames,'mean',sep='.')     

    pData.mean <- cbind(statpData,statistic='mean')
    rownames(pData.mean) <- colnames(exprs.mean)
    
    ## computing variances
    # we have to ensure sufficient replicates!
      which.enough.replicates <- which(statpData$.replicate>=nReplicatesVar)
      tmp.pData <- statpData[which.enough.replicates,,drop=FALSE]
      if (nrow(tmp.pData)>0){
        exprs.var <- do.call('cbind',
            lapply(igroups[which.enough.replicates],function(group) 
                apply(exprs(e)[,rownames(group)],1,var)))
        colnames(exprs.var) <- paste(irootnames[which.enough.replicates],'var',sep='.')     
            
        pData.var <- cbind(tmp.pData,statistic='var')
        rownames(pData.var) <- colnames(exprs.var)
      } else {
        exprs.var <- c()
        pData.var <- c()
      }

    ## computing diff to ref 
    reference.which <-  which(pData.mean$.reference)
    exprs.diffref <- exprs.mean[,-reference.which,drop=FALSE]- exprs.mean[,reference.which]
    colnames(exprs.diffref) <- paste(irootnames[-reference.which],'diffref',sep='.')
 
    pData.diffref <- cbind(statpData[!statpData$.reference,,drop=FALSE],statistic='diffref')
    rownames(pData.diffref) <- colnames(exprs.diffref)
  
        
    ## computing pooled standard deviations -- we have to check for replicates     
    # first ensure there are enough replicates for reference
  ### <TODO> add sign for pooledSD

    ref.replicates <- statpData[reference.which,'.replicates']
    if (ref.replicates<nReplicatesVar){
      warnings(paste('only ',ref.replicates,' replicates in the group ',iwithinname, ' to compute pooled standard errors( ',nReplicatesVar,'required)',sep=''))
      exprs.pooledSD <- c()
      pData.pooledSD <- c()
      exprs.spooledSD <- c()
      pData.spooledSD <- c()
      
    } else {
      # then select only groups with enough replicates
      tmp.pData <-  pData.var[pData.var$.replicates>=nReplicatesVar & !pData.var$.reference,,drop=FALSE]
      if (nrow(tmp.pData)>0){
        # ok, here there are some groups with enough replicates to compute pooled sd

        # formula:       # pooled variance =
        #   sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
 
        var1 <- exprs.var[,rownames(pData.var[pData.var$.reference,])]
        n1 <- pData.var[pData.var$.reference,'.replicates']
 
        exprs.pooledSD <- sapply(rownames(tmp.pData), 
          function(varvar){
            n2 <- tmp.pData[varvar,'.replicates']
            var2 <- exprs.var[,varvar]
            out <- sqrt(((n1-1)*var1 + (n2-1)*var2)/(n1+n2-2))
            return(out)
          })

        # then look at corresponding diff to reference and prepare signed SD for plot
        corresponding.diffref <- paste(substr(rownames(tmp.pData),1,nchar(rownames(tmp.pData))-4),'diffref',sep='.')
        exprs.spooledSD <- sapply(seq(along=rownames(tmp.pData)),
          function(ivarvar){
              varvar <- rownames(tmp.pData)[ivarvar]
              out <- exprs.pooledSD[,varvar]
#              diffname <- paste(tmp.pData[varvar,'.withingroups'],'diffref',sep='.')
#             diffname <- paste(irootnames[ivarvar],'diffref',sep='.')
            diffname <- corresponding.diffref[ivarvar]
              diffref <- exprs.diffref[,diffname]
              out[diffref<0] <- -out[diffref<0]
              return(out)
          })
      # prepare correct pData        
        cnames <- do.call('paste',
          c(as.list(tmp.pData[,c(within,reference$var,across),drop=FALSE]),'pooledSD',sep='.'))        
        colnames(exprs.pooledSD) <- cnames
        pData.pooledSD <- tmp.pData
        pData.pooledSD$statistic <- 'pooledSD'
        rownames(pData.pooledSD) <- colnames(exprs.pooledSD)
       

        cnames <- do.call('paste',
          c(as.list(tmp.pData[,c(within,reference$var,across),drop=FALSE]),'spooledSD',sep='.'))        
        colnames(exprs.spooledSD) <- cnames
        pData.spooledSD <- tmp.pData
        pData.spooledSD$statistic <- 'signedpooledSD'
        rownames(pData.spooledSD) <- colnames(exprs.spooledSD)
        
        
      }
      else {
        exprs.pooledSD <- c()
        pData.pooledSD <- c()
        exprs.spooledSD <- c()
        pData.spooledSD <- c()

      }
    }
    # we have done computations for this within / merge them to output object
    exprs.new  <- cbind(exprs.new,exprs.mean,exprs.var,exprs.diffref,exprs.pooledSD,exprs.spooledSD)
    pData.new  <- rbind(pData.new,pData.mean,pData.var,pData.diffref,pData.pooledSD,pData.spooledSD)
  }
  phenoData <- new('AnnotatedDataFrame',data=pData.new)
  out <- new('ExpressionSetWithComputation',
    exprs=as.matrix(exprs.new),
    phenoData=phenoData,
    featureData=featureData(e),
    experimentData=experimentData(e),
    annotation=annotation(e))
  invisible(out)  
}



# compute(e,reference=reference,within=within)

