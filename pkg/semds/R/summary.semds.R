summary.semds <- function(object, ...)
  {
    ## x ... object of class "semds"
    if (!is.null(object$thetatab)) {
      cat("\nTheta Parameters:\n")
      ttab <- cbind(round(object$thetatab[,1:2], 4), round(object$thetatab[,3:4], 3))
      print(ttab)
      cat("\n")
    } else {
      warning("No inference for saturated models.")
    }
  }
