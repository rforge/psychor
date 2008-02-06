`summary.smacofID` <-
function(object, ...)
{
  cat("\n")
  cat("Configurations:\n")
  print(lapply(object$conf, round, 4))
  #print(round(object$conf,4))
  cat("\nConfiguration dissimilarities: \n")
  print(lapply(object$confdiss, round, 4))
}


