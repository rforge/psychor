#active set methods for different solver


activeSet <- function(x, isomat, mySolver = lsSolver, ups = 1e-12, check = FALSE, ...) 
{
  
  
  a <- isomat
  extra <- list(...)
  n <- length(x)
  xold <- x
  ax <- aTx(a, xold)
  ia <- is.active(ax, ups=ups)

  repeat {
    if (length(ia)==0) { 
      aia<-NULL
    } else { 
      aia<-a[ia,]  
    }
    yl <- mySolver(xold, aia, extra)
    
    y <- yl$y
    lbd <- yl$lbd
    fy <- yl$f
    gy <- yl$gy
    ay <- aTx(a,y)
    iy <- which.min(ay)
    my <- ay[iy]
    if (length(lbd)==0) {
      ml <- Inf
    } else {
      il <- which.min(lbd)
      ml <- lbd[il]
    }
    
    if (is.pos(my,ups)) {
      if (is.pos(ml,ups)) break()
      xnew <- y
      ax <- ay
      ia <- ia[-il]
    } else {
      k <- which((ax>0)&(ay<0))
      rat <- -ay[k]/(ax[k]-ay[k])
      ir <- which.max(rat)
      alw <- rat[ir]
      xnew <- y+alw*(xold-y)
      ax <- aTx(a,xnew)
      ia <- sort(c(ia,k[ir]))
    }
    xold <- xnew
  }

  lup <- rep(0, length(ay)) 
  lup[ia] <- lbd
  hl <- taTx(a, lup, n)
  if (check) {
    ck <- checkSol(y, gy, a, ay, hl, lup, ups)
  } else {
    ck <- NULL
  }

  result <- list(x = y, lbd = lup, fval = fy, ay = ay, hl = hl, gradient = gy, isocheck = ck, call = match.call())
  class(result) <- "activeset"
  result
}