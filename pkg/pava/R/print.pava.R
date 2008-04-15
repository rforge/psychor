`print.pava` <-
function(x, ...)
{
  #x... object of class pava
  cat("Fitted values: \n")
  print(x$yfit)
  cat("\n")
}

