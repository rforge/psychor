plot.gpava <- function (x, main = "PAVA Plot", xlab = "Predictor", ylab = "Response", ...)
{
#x ... object of class pava

  o <- order(x$x)
  xval <- x$x[o]
  yval <- x$yfit[o]
  xcum <-  c(xval[1] - mean(diff(xval)), xval) 
  jumps <- ((1:length(yval))[!duplicated(yval)]-1)[-1]   #jumps of fitted step function
 
  if (is.list(x$y)) {
    plot(xcum, c(NA, (mapply(x$solver, x$y, x$w)[o])), xlab = xlab, ylab = ylab, main = main, ...) 
  } else {
    plot(xcum, c(NA, x$y[o]), xlab = xlab, ylab = ylab, main = main, ...)   
  }  
  lines(xval, yval, col = "lightblue", lwd = 1, type = "S")
  points(xval[jumps], yval[jumps], col = "lightblue", pch = 13)
  grid()
}
