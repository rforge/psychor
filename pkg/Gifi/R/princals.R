princals <- function(data, ndim = 2, level = "ordinal", 
                     active = TRUE, eps = 1e-6, itermax = 1000, verbose = 0)
## NLPCA wrapper function for Gifi 
{
  
  name <- deparse(substitute(data))  
  call <- match.call()
  if (ndim > ncol(data)) stop("Number of dimensions can't be larger than the number of variables!")
  
  ## ---------------------------------- Gifi call -------------------------------------------
  fit <- gifiEngine(data = data, ndim = ndim, rank = 1, level = level, active = active, eps = eps, 
                itermax = itermax, verbose = verbose, call = call)        ## rank-1 call
  
  
  ## --------------------------------- PCA-like outputs ------------------------------------
  scores1 <- fit$scoremat[,,1]                               ## matrix with category scores (first dimension)
  R1 <- cor(scores1)                                         ## new correlation matrix
  eigenvalsall <- eigen(R1)$values                           ## all eigenvalues
  eigenvals <- eigenvalsall[1:ndim]
  varexp <- eigenvals/(sum(eigenvalsall))*100                ## amount of explained variance
  
  loads <- fit$loadings                                      ## loadings
  loadsmat <- as.matrix(do.call(rbind, loads))
  rownames(loadsmat) <- colnames(data)
  loadstand <- apply(loadsmat, 2, function(ll) ll/sqrt(sum(ll^2)))  ## standardized loadings
  colnames(loadstand) <- paste0("Comp.", 1:ncol(loadstand))
  
  pcscores <- scores1 %*% loadstand                          ## principal component scores
    
  ## colSums(loadsmat^2)*100   ## same as through eigenvalues
  
  result <- list(datname = name, eigenvalues = eigenvals, vaf = varexp, loadings = loadstand, pcscores = pcscores, 
                 catscores = fit$catscores, scoremat = fit$scoremat, objscores = fit$objscores, 
                 cat.centroids = fit$cat.centroids, ind.mat = fit$ind.mat,  
                 discrim = fit$discrim, ndim = fit$ndim, niter = fit$niter, level = fit$level, 
                 loss = fit$loss, active = fit$active, dframe = fit$dframe, call = match.call())
  class(result) <- c("princals", "gifi")
  return(result)
}