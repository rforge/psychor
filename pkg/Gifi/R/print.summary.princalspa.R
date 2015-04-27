print.summary.princalspa <- function(x, ...)
{
  cat("\n")
  cat("Results parallel analysis:\n")
  print(round(x$partab, 3))
  cat("\n")
}