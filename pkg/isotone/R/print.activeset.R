print.activeset <- function(x, ...)
{
  cat("\nValue target function: \n")
  print(x$fval)
  cat("\nFitted values: \n")
  print(x$x)
  cat("\n")
}