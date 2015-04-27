print.princalsrot <- function(x, ...) {
  
  cat("Call: ")
  print(x$call)
  
  cat("\nLoadings after", toupper(x$method), "rotation:\n\n")
  print(round(x$loadings, 3))
  cat("\n")
}