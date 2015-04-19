print.summary.princals <- function(x, ...)
{
  cat("\nImportance of Components:\n")
  print(x$vartab)
  cat("\n")
}