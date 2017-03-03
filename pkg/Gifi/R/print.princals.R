print.princals <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Loss value:", round(x$f, 3),"\n")
  cat("Number of iterations:",x$ntel,"\n")
  cat("\n")
  cat("Importance of Components:\n")
  nvars <- ncol(x$loadings)
  ndim <- nrow(x$lambda)
  df <- as.data.frame(rbind(round(x$evals[1:ndim], 4), 
                            round((x$evals/nvars*100)[1:ndim], 4), 
                            round(cumsum((x$evals/nvars)[1:ndim])*100, 2)))
  colnames(df) <- paste0("Comp", 1:ncol(df))
  rownames(df) <- c("Eigenvalues", "VAF", "Cumulative VAF")
  print(df)
  cat("\n")
}
