
addQuantilesColors <- function(e,ngroups=3)
{
#D ngroups=5
### will add quantiles levels in ExpressionSet featureData slot
  stopifnot(is(e,'ExpressionSet'))
  means <- apply(exprs(e),1,mean,na.rm=TRUE)
  quantiles <- quantile(means,probs=seq(0,1,length=ngroups+1))
  colorsVector <- cut(means,quantiles,include.lowest=TRUE,labels=FALSE)
  if ('colorsQuantilesVector' %in% colnames(fData(e))) 
    fData(e)[,'colorsQuantilesVector'] <- colorsVector else 
    fData(e) <- data.frame(cbind(as.matrix(fData(e)),colorsQuantilesVector=colorsVector))
    invisible(e)    
  }
