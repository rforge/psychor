# update function for transformation of external variables
updext <- function(x,w,external,extvars,constraint){
  m <- ncol(external)
  # reconstruct weigh matrix C
  if (constraint == "linear"){
    C <- solve(crossprod(external,w%*%external),crossprod(external,w%*%x))
  } else if (constraint == "diagonal") {
    C <- diag(colSums(external*(w%*%x))/colSums(external*(w%*%external)))
  }
  # For updating external[,s] we need a value larger than the largest eigenvalue
  # of V kronecker CC'.  
  svdC <- svd(C)
  v2 <- 2*max(diag(w))* svdC$d[1]^2                     # 2*max(diag(w)) is an upperbound of the largest eigenvalue of w
  z <- external - (1/v2)*(w %*% (external %*% C - x)) %*% t(C) # Compute the unconstrained majorization update
  external.old <- external
  iord.prim <- list()
  for (s in 1:ncol(external)){
    tt <- transform(z[,s], extvars[[s]], normq = 0)     # Compute update for external variable s
    external[,s] <- tt$res*(n/sum(tt$res^2))^.5         # Make the external variable of length n
    iord.prim[[s]] <- tt$iord.prim                      # Retain the ordening if primary approach to ties
  }
  return(list(external = external, iord.prim = iord.prim))
}
