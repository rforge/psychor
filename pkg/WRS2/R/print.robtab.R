print.robtab <-
function(x,...)
{
  cat("Call:\n")
  print(x$call)
  
  cat("\nParameter table: \n")
  print(x$partable)
  cat("\n")
}
