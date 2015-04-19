morals <- function(data, resp = 1, level = "ordinal", active = TRUE, 
                     eps = 1e-6, itermax = 1000, verbose = 0)
## MORALS wrapper function for homals 
{
  name <- deparse(substitute(data))  
  cols <- seq(1, ncol(data))
  
  fit <- gifi(data = data, ndim = 1, sets = list(c(resp), cols[-resp]),
              level = level, active = active, eps = eps, 
              itermax = itermax, verbose = verbose)          ## dim-1 call
  

  result <- list(datname = name, catscores = fit$catscores, scoremat = fit$scoremat, objscores = fit$objscores, 
                 cat.centroids = fit$cat.centroids, ind.mat = fit$ind.mat, loadings = fit$loadings, 
                 low.rank = fit$low.rank, discrim = fit$discrim, ndim = fit$ndim, niter = fit$niter, level = fit$level, 
                 eigenvalues = fit$eigenvalues, loss = fit$loss, rank.vec = fit$rank.vec, active = fit$active, 
                 dframe = fit$dframe, call = match.call())
  
  class(result) <- c("morals", "gifi")
  return(result)
}