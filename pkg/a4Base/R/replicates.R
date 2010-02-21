`replicates` <-
function(x)
{
  x <- as.character(x)
  x[is.na(x)] <- ''
  replicate <- vector("numeric", len=length(x))
  replicate[order(x)] <- unlist(sapply(rle(sort(x))$lengths, seq_len))
  return(replicate)
}

