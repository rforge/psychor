plot.pava <- function (x, main = "PAVA Plot", xlab = "Predictor", ylab = "Response", ...)
{
#x ... object of class pava


# plot.type = c("single", "row.wise", "col.wise")
# par.fit = list(col = "red", cex = 1.5, pch = 13, lwd = 1.5)
# mar = if (both) 0.1 + c(3.5, 2.5, 1, 1) else par("mar"),
# mgp = if (both) c(1.6, 0.7, 0) else par("mgp"),
# grid = length(x$x) < 12

  o <- order(x$x)
  xval <- x$x[o]
  xcum <-  c(xval[1] - mean(diff(xval)), xval) 
  jumps <- ((1:length(x$yfit))[!duplicated(x$yfit)]-1)[-1]   #jumps of fitted step function
 
  plot(xcum, c(NA, x$y), xlab = xlab, ylab = ylab, main = main, ...)   
  lines(xval, x$yfit, col = "lightblue", lwd = 1, type = "S")
  points(xval[jumps], x$yfit[jumps], col = "lightblue", pch = 13)
  grid()
}
