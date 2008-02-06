`print.smacofID` <-
function(x,...)
{
  cat("\n")
  cat("Model:",x$model,"\n")
  cat("Number of objects:",x$nobj,"\n")
  cat("\nStress:",x$stress.m,"\n")
  if(x$stress.m != x$stress.nm) cat("Nonmetric stress:",x$stress.nm,"\n")
  if(x$stress.m != x$stress.uc) cat("Unconstrained stress",x$stress.uc,"\n")
  if(x$stress.m != x$stress.c) cat("Constrained stress",x$stress.co,"\n")
   
  cat("Number of iterations:",x$niter,"\n")
  cat("\n")
}


