print.princals <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Loss value:", round(x$f, 3),"\n")
  cat("Number of iterations:",x$ntel,"\n")
  cat("\n")
  cat("Eigenvalues:", round(x$evals[1:ncol(x$loadings)], 3), "\n")
  cat("\n")
}
