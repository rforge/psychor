summary.princalspa <- function(object, ...) {
  
  partab <- data.frame(PCeigen = object$pceigen, Simeigen = object$simeigen, sd = object$simsd, 
                       lowerCI = object$simci[1,], upperCI = object$simci[2,])
  res <- list(partab = partab)
  class(res) <- "summary.princalspa"
  return(res)
}