print.aspect <- function(x, ...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("\nCorrelation matrix of the scaled data:\n")
  print(x$cormat)
  cat("\n")
  invisible(x)
}
