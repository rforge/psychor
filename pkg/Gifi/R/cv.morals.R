cv.morals <- function(object, verbose = FALSE, ...) {
  ## object of class morals
  if ((class(object)[1]) != "morals") stop("CV implemented for Morals only!")
  
  n <- length(object$yhat)
  ypreds <- ypreds <- numeric(n)
  for (i in 1:n) {
    if (verbose) cat("Left-out observation: ", i, "\n")
    dat <- object$data[-i,]
    morcall <- object$call
    morcall$x <- dat[,-1]
    morcall$y <- dat[, 1]             
    morres <- eval(morcall)          ## re-fit morals
    xleft <- object$xhat[i,]         ## left-out observation
    class(xleft) <- "numeric"
    qxy <- lsRC(morres$xhat, morres$yhat)$solution   ## coefficient
    ypreds[i] <- xleft %*% qxy                       ## predict left-out observation
  }
  yobs <- object$yhat
  ei <- ypreds - yobs
  sqrt(mean((ei)^2))    ## RMSE prediction
} 


