print.summary.princals <- function(x, ...)
{
  cat("\nLoadings: \n")
  print(round(x$loadings, 3))
  cat("\nVariance accounted for:\n")
  print(round(x$vartab, 2))
  cat("\n")
}