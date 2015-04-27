print.princalspa <- function(x, ...) {
  
  cat("Call: ")
  print(x$call)
  
  cat("\nThe suggested number of components is ", x$ncomp, ".\n\n", sep = "")
  
}