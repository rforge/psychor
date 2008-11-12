print.activeset <- function(x, ...)
{
  cat("\nValue target function: \n")
  print(x$func.vals)
  cat("\nFitted values: \n")
  print(x$x)
  cat("\n")
}