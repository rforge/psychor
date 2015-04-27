print.summary.princalsvss <- function(x, ...)
{
  cat("\n")
  cat("Results VSS and MAP analysis:\n")
  print(round(x$vsstab, 3))
  cat("\n")
}