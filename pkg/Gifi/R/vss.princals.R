vss <- function(x, ndim = NULL) UseMethod("vss")

## very simple structure for princals objects
vss.princals <- function(x, ndim = NULL) {

  if (is.null(ndim)) {
    ndim <- x$ndim
  } else {
    call <- x$call
    call$ndim <- ndim                           ## assign the updated dimensions
    x <- eval(call, parent.frame())             ## re-fit princals
  }
  
  ## ------- perform VSS
  resvss <- psych::vss(x$scoremat[,,1], n = ndim, rotate = "none", fm = "pc", plot = FALSE)
    
  res <- list(VSS1 = resvss$cfit.1, VSS2 = resvss$cfit.2, MAP = resvss$map, sqresid = resvss$vss.stats$sqresid, 
                      fit = resvss$vss.stats$fit, call = match.call())
  class(res) <- "princalsvss"
  res
}