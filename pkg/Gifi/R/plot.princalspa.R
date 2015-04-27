plot.princalspa <- function(x, main = "Parallel Analysis", xlab = "Number of Components", 
                             ylab = "Eigenvalues", leg.pos = "topright", leg.lab = c("Fitted", "Simulated"), ylim, ...) {
  nc <- length(x$pceigen)
  if (missing(ylim)) ylim <- c(0, max(x$pceigen))
  
  plot(1:nc, x$pceigen, type = "b", pch = 19, ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)
  points(1:nc, x$simeigen, type = "b", lty = 2)
  legend(leg.pos, legend = leg.lab, lty = c(1,2))
}