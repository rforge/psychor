plot.princals <- function(x, plot.type = "biplot", col.scores = "gray", col.loadings = "black", arrows = TRUE, cex.scores = 0.8, 
                          cex.loadings = 0.8, labels.scores = FALSE, labels.loadings = TRUE, var.subset = "all", plot.dim = 1, 
                          main, xlab, ylab, xlim, ylim, ...)
{

## S3 plot method for objects of class "princals"
## Produces various 2D-plots
## plot.dim ... vector of length 2 with dimensions to be plotted against
## plot.type ... type of plot to be drawn: "loadplot", "screeplot", "biplot"

if ((x$ndim == 1) && (plot.type != "transplot")) stop("No biplot/loadings plot can be drawn for ndim = 1!")
nvar <- dim(x$dframe)[2]
match.arg(plot.type, c("biplot", "loadplot", "screeplot", "transplot"))


#----------------------------------loadplot-------------------------------------
if (plot.type == "loadplot") {
  xycoor <- x$loadings[, 1:2]
  if (missing(xlim)) {                            ## x limits
     xlim.min <- min(xycoor[,1],0)
     xlim.max <- max(xycoor[,1],0)
     xlim <- c(xlim.min, xlim.max)*1.2
  }
  if (missing(ylim)) {                            ## y limits
    ylim.min <- min(xycoor[,2],0)
    ylim.max <- max(xycoor[,2],0)
    ylim <- c(ylim.min, ylim.max)*1.2
  }
  if (missing(xlab)) xlab <- "Dimension 1"       ## labels
  if (missing(ylab)) ylab <- "Dimension 2"
  if (missing(main)) main <- "Loadings Plot"
  
  plot(xycoor, type = "p", pch = ".", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main, 
       cex = cex.loadings, col = col.loadings, ...)
  abline(h = 0, v = 0, col = "gray", lty = 2)
  for (i in 1:nvar) arrows(0, 0, xycoor[i,1],xycoor[i,2], length = 0.08)   #lines(rbind(xycoor[i,],c(0,0)))
  posvec <- apply(xycoor, 1, sign)[2,] + 2      
  text(xycoor, labels = rownames(xycoor), pos = posvec, cex = cex.loadings)
    
}
#-------------------------------- end loadplot ---------------------------------

## --------------------------------- biplot ------------------------------------
if (plot.type == "biplot") {
  xycoorL <- x$loadings[, 1:2]
  xycoorS <- x$pcscores[, 1:2]
  if (missing(xlim)) {                            ## x limits
    xlim.min <- min(c(xycoorL[,1], xycoorS[,1]), 0)
    xlim.max <- max(c(xycoorL[,1], xycoorS[,1]), 0)
    xlim <- c(xlim.min, xlim.max)*1.05
  }
  if (missing(ylim)) {                            ## y limits
    ylim.min <- min(c(xycoorL[,2], xycoorS[,2]), 0)
    ylim.max <- max(c(xycoorL[,2], xycoorS[,2]), 0)
    ylim <- c(ylim.min, ylim.max)*1.05
  }
  if (missing(xlab)) xlab <- "Dimension 1"       ## labels
  if (missing(ylab)) ylab <- "Dimension 2"
  if (missing(main)) main <- "Biplot"
  
  if (labels.scores) {
    plot(xycoorS, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)
    text(xycoorS, labels = rownames(xycoorS), col = col.scores, cex = cex.scores)
  } else {
    plot(xycoorS, type = "p", pch = 19, cex = cex.scores, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main, 
         col = col.scores, ...)
  }
  points(xycoorL, pch = ".")
  abline(h = 0, v = 0, col = "gray", lty = 2)
  for (i in 1:nvar) arrows(0, 0, xycoorL[i,1],xycoorL[i,2], length = 0.08)   #lines(rbind(xycoor[i,],c(0,0)))
  posvec <- apply(xycoorL, 1, sign)[2,] + 2      
  text(xycoorL, labels = rownames(xycoorL), pos = posvec, cex = cex.loadings) 
}
#------------------------------------trfplot------------------------------------
#draws transformation plots
## FIXME: use lattice

if (plot.type == "transplot") {
 
    if (missing(main)) main <- colnames(x$dframe)
    if (missing(xlab)) xlab <- "Observed"
    if (missing(ylab)) ylab <- "Transformed"
    
    if (var.subset[1] == "all") var.subset <- colnames(x$dframe)       ## extract variables and scores to be plotted
    nvars <- length(var.subset)                                 ## number of variables to be plotted
    plotvars <- x$dframe[,var.subset]
    plotcats <- lapply(plotvars, levels)
    ploty <- x$catscores[var.subset]
    
    
    ## set up number of vertical and horizontal panels  
    npanv <- ceiling(sqrt(nvars)) 
    npanh <- floor(sqrt(nvars))
    if (npanv * npanh < nvars) npanv <- npanv + 1
    
    op <- par(mfrow = c(npanv, npanh))
    for (i in 1:nvars) {
      
      w <- options("warn")                               ## x-values
      options(warn = -1)
      xvals <- as.numeric(plotcats[[i]])
      if (any(is.na(xvals))) xvals <- 1:length(plotcats[[i]])
      options(w)
      
      yvals <- ploty[[i]][, plot.dim]                    ## y-values
      
      if (length(xvals) > 20) pch <- "." else pch <- 19
      plot(xvals, yvals, pch = pch, type = "b", xlab = xlab, ylab = ylab, main = main[i], xaxt = "n", ...)
      axis(1, at = xvals, labels = plotcats[[i]])
    }
    par(op)

}    
    # #----------------------------------end trfplot----------------------------------


#---------------------------------- screeplot ----------------------------------
if (plot.type == "screeplot") {

   if (missing(main)) main <- "Scree Plot"
   if (missing(xlab)) xlab <- "Number of Dimensions"
   if (missing(ylab)) ylab <- "Eigenvalues"
  
   nd <- length(x$eigenvalues)
   plot(1:nd, x$eigenvalues, type = "b", xlab = xlab, ylab = ylab, main = main, xaxt = "n", pch = 19)
   axis(1, at = 1:nd, labels = 1:nd)
}  
#-------------------------------- end screeplot --------------------------------
 
}

