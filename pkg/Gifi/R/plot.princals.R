plot.princals <- function(x, plot.dim = c(1,2), plot.type = "loadplot", var.subset, main, type, xlab, ylab, 
         xlim, ylim, identify = FALSE, ...)
{
## FIXME: biplot
## S3 plot method for objects of class "princals"
## Produces various 2D-plots
## plot.dim ... vector of length 2 with dimensions to be plotted against
## plot.type ... type of plot to be drawn: "loadplot", "screeplot", "trfplot"

options(locatorBell = FALSE)
if (x$ndim == 1) stop("No plots can be drawn for ndim = 1!")
if (length(plot.dim) !=  2) stop("plot.dim must be of length 2!")
if ((plot.type != "trfplot") && (plot.type != "screeplot")) {      #plot.dim are ignored for trfplot
  pd1 <- plot.dim[1]
  pd2 <- plot.dim[2]
  if (pd2 > x$ndim) stop("Only",x$ndim,"dimensions were extracted!")
  if (missing(xlab)) xlab <- paste("Dimension",pd1)
  if (missing(ylab)) ylab <- paste("Dimension",pd2)
}

nvar <- dim(x$dframe)[2]
if (missing(var.subset)) var.subset <- 1:nvar
   

#----------------------------------loadplot-------------------------------------
if (plot.type == "loadplot") {
  xycoor <- x$loadings[, c(pd1,pd2)]
  if (missing(main)) main1 <- "Loadings Plot" else main1 <- main

  xlim.min <- min(xycoor[,1],0)
  xlim.max <- max(xycoor[,1],0)
  ylim.min <- min(xycoor[,2],0)
  ylim.max <- max(xycoor[,2],0)
  if (missing(xlim)) xlim <- c(xlim.min,xlim.max)*1.2
  if (missing(ylim)) ylim <- c(ylim.min,ylim.max)*1.2

  plot(xycoor,type = "p", pch = 19, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, main = main1, cex = 0.5, ...)
  abline(h = 0, v = 0, col = "gray", lty = 2)
  for (i in 1:nvar) arrows(0, 0, xycoor[i,1],xycoor[i,2], length = 0.08)   #lines(rbind(xycoor[i,],c(0,0)))
  abline(h = 0, col = "lightgray", lty = 2)
  abline(v = 0, col = "lightgray", lty = 2)
  
  if (identify) {
    identify(xycoor, labels = rownames(xycoor), cex = 0.8)
  } else {
    posvec <- apply(xycoor, 1, sign)[2,] + 2      
    text(xycoor, labels = rownames(xycoor), pos = posvec, cex = 0.8)
  }  
}
#-------------------------------- end loadplot ---------------------------------




#------------------------------------trfplot------------------------------------
#draws transformation plots
## FIXME: use lattice

# if (plot.type == "trfplot") {
# 
#    if (missing(type)) type <- "l"
#    if (missing(xlab)) xlab <- "original scale"
#    if (missing(ylab)) ylab <- "transformed scale"
# 
#    for (i in var.subset) {
#      
#      if (missing(main)) main1 <- paste("Transformation plot for", colnames(x$dframe[i])) else main1 <- main
#      if (missing(ylim)) ylim1 <- range(x$low.rank[[i]]) else ylim1 <- ylim
#      
#      p <- dim(x$low.rank[[i]])[2]                           #number of dimensions
#      vlev <- rownames(x$low.rank[[1]])
#        
#      op <- par("ask" = TRUE)                     #first dimensions
#      matplot(x$low.rank[[i]], type = type, main = main1, ylim = ylim1, xlab = xlab, 
#             ylab = ylab, xaxt = "n", pch = 19)
#      if (p != 1) legend(leg.pos,paste("Solution",1:p),col = 1:p, lty = 1:p)
#      axis(1, at = 1:dim(x$low.rank[[i]])[1], labels = rownames(x$low.rank[[i]]))
#      
#    }
# }
# #----------------------------------end trfplot----------------------------------


#---------------------------------- screeplot ----------------------------------
if (plot.type == "screeplot") {

   if (missing(main)) main <- "Scree Plot"
   if (missing(xlab)) xlab <- "Dimension"
   if (missing(ylab)) ylab <- "Eigenvalue"
  
   nd <- length(x$eigenvalues)
   plot(1:nd, x$eigenvalues, type = "b", xlab = xlab, ylab = ylab, main = main, xaxt = "n", pch = 19)
   axis(1, at = 1:nd, labels = 1:nd)
}  
#-------------------------------- end screeplot --------------------------------
 
}

