svm_mdsplot <- function(mds_object, svm_object, class, legend1 = TRUE, legend2 = TRUE, inset = c(-0.2, 0.5)) {
  X <- mds_object$conf
  nc <- svm_object$nclasses
  
  xlim <- range(X[, 1])*1.1
  ylim <- range(X[, 2])*1.1
  minmax <- range(c(X))*2.2      
  seqs <- seq(minmax[1], minmax[2], 0.01)
  
  dnames <- attr(svm_object$terms, "term.labels")
  xgrid <- expand.grid(list(seqs, seqs))
  colnames(xgrid) <- dnames
  lev <- predict(svm_object, newdata = xgrid)
  levmat <- matrix(as.numeric(lev), sqrt(length(lev)))
  op <- par(pty = "s", mar = c(5.1, 4.1, 4.1, 9.1), xpd = TRUE)
  image(seqs, seqs, levmat, xlab = "D1", ylab = "D2", col = rainbow_hcl(nc), asp = 1,
        xlim = xlim, ylim = ylim, main = "SVM Facets")
  points(X, cex = 0.7, pch = (1:nc)[as.numeric(class)])      
  text(X, labels = rownames(X), pos = 3, cex = 0.7)
  if (legend1) legend("topright", inset = c(inset[1], 0), legend = svm_object$levels, pch = 15, title = "Facets", col = rainbow_hcl(nc), cex = 0.7)
  if (legend2) legend("topright", inset = c(inset[1], inset[2]), legend = svm_object$levels, title = "Classes", pch = 1:nc, cex = 0.7)
  par(op)
}
