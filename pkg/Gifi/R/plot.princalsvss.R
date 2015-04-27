plot.princalsvss <- function(x, main = "Very Simple Structure", xlab = "Number of Components", 
                             ylab = "VSS", ylim = c(0,1), leg.pos = "bottomleft", leg.lab = c("VSS1", "VSS2"), ...) {
  nc <- length(x$VSS1)
  plot(1:nc, x$VSS1, type = "b", pch = 19, ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)
  points(2:nc, x$VSS2[-1], type = "b", lty = 2)
  legend(leg.pos, legend = leg.lab, lty = c(1,2))
}