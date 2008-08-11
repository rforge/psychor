print.aspect <- function(x, ...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("\nLoss:",x$loss)
  cat("\n\nCorrelation matrix of the transformed data:\n")
  print(x$cormat)
  cat("\n")
  invisible(x)
}
