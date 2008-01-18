`print.smacofB` <-
function(x,...)
{
  cat("\n")
  cat("Model:",x$model,"\n")
  cat("Number of objects:",x$nobj,"\n")
  cat("\nMetric stress:",x$stress.m,"\n")
  if(x$stress.m != x$stress.nm) cat("Nonmetric stress:",x$stress.nm,"\n")
  cat("Number of iterations:",x$niter,"\n")
  cat("\n")
}

