parallel <- function(x, N = 1000, method = "permutation") UseMethod("parallel")
  
## parallel analysis for princals objects
parallel.princals <- function(x, N = 1000, method = "permutation") {
    
  match.arg(method, c("random", "permutation"))  
  
  if (x$ndim != ncol(x$scoremat[,,1])) {        ## in case the full model was not fitted
    ndim <- ncol(x$scoremat[,,1])
    call <- x$call
    call$ndim <- ndim                           ## assign the updated dimensions
    x <- eval(call, parent.frame())             ## re-fit princals
  }
  
  X <- x$scoremat[,,1]
  n <- nrow(X)
  p <- ncol(X)
     
  evmat <- matrix(NA, N, p)
  if (method == "random") {                    ## random (normal) data matrices
    for (i in 1:N) 
    {
      R <- cor(matrix(rnorm(n * p), n, p))
      evals <- eigen(R)$values
      evmat[i,] <- evals
    }
  }
  if (method == "permutation") {
    for (i in 1:N) {
       Xi <- apply(X, 2, function(cols) {
               cols[sample(1:n)]                   ## permute within columns
       })
       R <- cor(Xi)
       evals <- eigen(R)$values
       evmat[i,] <- evals
    }
  }
    
  evalsSim <- colMeans(evmat)                        ## mean eigenvalues 
  sdSim <- apply(evmat, 2, sd)                       ## sd eigenvalues
  ciSim <- apply(evmat, 2, quantile, probs = c(0.025, 0.975))  ## CI eigenvalues
  colnames(ciSim) <- paste0("Comp.", 1:p)
  
  ncomp <- sum(x$eigenvalues > evalsSim)
  
  res <- list(pceigen = x$eigenvalues, simeigen = evalsSim, simsd = sdSim, simci = ciSim, ncomp = ncomp, call = match.call())
  class(res) <- "princalspa"
  res
}
  
