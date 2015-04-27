Gifi <-
function(data, ndim = 2, rank = ndim, level = "nominal", sets = 0, active = TRUE, 
eps = 1e-6, itermax = 1000, verbose = 0)
{
#data ... data frame
#sets ...  list of vectors of set indices 
#level ... which measurement level (either single string or vector
#ndim ... number of dimensions
#active ... which variables are active (single TRUE means all)
#rank ... which category quantification ranks (default all ndim)
#eps ... iteration precision eigenvalues (default 1e-6)

name <- deparse(substitute(data))  
call <- match.call()
fit <- gifiEngine(data = data, ndim = ndim, rank = rank, level = level, sets = sets, active = active, eps = eps, 
               itermax = itermax, verbose = verbose, call = call)        

result <- list(datname = name, catscores = fit$catscores, scoremat = fit$scoremat, objscores = fit$objscores, 
               cat.centroids = fit$cat.centroids, ind.mat = fit$ind.mat, loadings = fit$loadings, 
               low.rank = fit$low.rank, discrim = fit$discrim, ndim = fit$ndim, niter = fit$niter, level = fit$level, 
               eigenvalues = fit$eigenvalues, loss = fit$loss, rank.vec = fit$rank.vec, active = fit$active, 
               dframe = fit$dframe, call = match.call())
class(result) <- "gifi"
result
}

