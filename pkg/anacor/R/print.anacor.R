`print.anacor` <-
function(x, ...)
{
# x ... object of class "anacor"
  cat("\nCA fit: \n")
  cat("\nTotal chi-square value:",round(x$chisq,3),"\n")
  cat("Sum of eigenvalues (total inertia):", round(sum(x$eigen.values), 3),"\n")
  cat("Eigenvalues (principal inertias):\n")
  cat(round(x$eigen.values, 3),"\n")
  cat("\n")
  if (!is.null(x$rmse[[1]])) cat("Benzecri RMSE rows: ", x$rmse[[1]],"\n")
  if (!is.null(x$rmse[[2]])) cat("Benzecri RMSE columns: ", x$rmse[[2]],"\n")
  
  cat("\nChi-square decomposition: \n")
  print(round(x$chisq.decomp,3))
  
  if (!is.null(x$stestmat)) {
   cat("\nz-test on singular values: \n")
   print(round(x$stestmat, 3))
  }
  cat("\n")
}

