rotation <- function(x, method = "varimax", normalize = TRUE, eps = 1e-5, m = 4) UseMethod("rotation")

## varimax/promax rotation for princals
rotation.princals <- function(x, method = "varimax", normalize = "TRUE", eps = 1e-5, m = 4) {
  
  ## -------------------------------------- rotation ----------------------------
  ndim <- x$ndim
  if (ndim == 1) stop("No rotation performed for ndim = 1.")
  match.arg(method, c("promax", "varimax"))
    
  p <- nrow(x$loadings)                                               ## number of variables
  rawload <- x$loadings %*% diag(sqrt(x$eigenvalues), ndim, ndim)    ## raw loadings (eigenvectors scaled by sqrt of eigenvalues)
  
  if (method == "varimax") {
    rotres <- stats::varimax(rawload, normalize = normalize, eps = eps)  ## varimax rotation               
    rotload <- unclass(rotres$loadings)                                  ## extract loadings
    rotload <- rotload %*% solve(diag(sqrt(x$eigenvalues), ndim, ndim))      ## scale back loadings
    rotscores <- x$scoremat[,,1] %*% t(MASS::ginv(rotload))
    rotmat <- rotres$rotmat
  }  
  
  if (method == "promax") {
    rotres <- stats::promax(rawload, m = m)                   ## promax rotation               
    rotload <- unclass(rotres$loadings)                                  ## extract loadings
    rotload <- rotload %*% solve(diag(sqrt(x$eigenvalues), ndim, ndim))      ## scale back loadings
    rotscores <- x$scoremat[,,1] %*% t(MASS::ginv(rotload))
    rotmat <- rotres$rotmat
  }
  
  colnames(rotload) <- colnames(rotscores) <- paste0("Comp.", 1:ndim)
  
  res <- list(loadings = rotload, pcscores = rotscores, rotmat = rotmat, method = method, call = match.call())
  class(res) <- "princalsrot"
  return(res)
} 
  
