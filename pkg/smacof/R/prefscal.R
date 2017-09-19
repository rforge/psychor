`prefscal` <- function(delta, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"),
                       conditionality = c("matrix", "row"), lambda = 0.5, omega = 0.1, 
                       circle = c("none", "row", "column"), weightmat = NULL, init = NULL, 
                       ties = c("primary","secondary", "tertiary"), verbose = FALSE, itmax = 1000, reg = 1e-6, 
                       eps = 1e-6, spline.degree = 2, spline.intKnots = 2)
  
  # init ... either a list of 2 matrices of dimension n \times p and m \times p with 
  # starting values. if NULL, svd is used.
{
  circle <- match.arg(circle, c("none", "row", "column"), several.ok = FALSE)
  type <- match.arg(type, c("ratio", "interval", "ordinal", "mspline"), several.ok = FALSE)
  ties <- match.arg(ties, c("primary","secondary", "tertiary"), several.ok = FALSE)
  conditionality <- match.arg(conditionality, c("matrix", "row"), several.ok = FALSE)
  diss <- delta
  rnames <- rownames(delta)
  if (is.data.frame(diss)) diss <- as.matrix(diss)
  checkdiss(diss)
  
  n <- dim(diss)[1]                       #number of individuals
  m <- dim(diss)[2]                       #number of objects
  p <- ndim
  
  if (is.null(weightmat)) {
    w <- matrix(1,n,m)                    #initialize weights (as 1)
  } else w <- weightmat
  sumw <- sum(w)
  
  delta <- ifelse(is.na(diss),0,diss)     #replace NA's by 0
  
  ## --- Prepare for optimal scaling
  trans <- type
  if (trans=="ratio"){
    trans <- "none"
  } else if (trans=="ordinal" & ties=="primary"){
    trans <- "ordinalp"
  } else if(trans=="ordinal" & ties=="secondary"){
    trans <- "ordinals"
  } else if(trans=="ordinal" & ties=="tertiary"){
    trans <- "ordinalt"
  } else if(trans=="spline"){
    trans <- "mspline"
  }
  disobj <- list()
  if (conditionality == "matrix"){
    disobj[[1]] <- transPrep(as.vector(delta), trans = trans, spline.intKnots = spline.intKnots, spline.degree = spline.degree)
    tt <- transform(as.vector(delta), disobj[[1]], w = as.vector(w), normq = sumw)
    dhat <- matrix(tt$res, n, m)  ## dhat update
  } else { ## conditionality == "row"
    dhat <- matrix(0, n, m)
    for (i in 1:n) {
      disobj[[i]] <- transPrep(delta[i, ], trans = trans, spline.intKnots = spline.intKnots, spline.degree = spline.degree)
      dhat[i, ] <- transform(delta[i, ], disobj[[i]], w = w[i, ], normq = m)$res  ## dhat update
    }
  }
  ## --- end optimal scaling prep
  
  
  itel <- 1
  #delta <- delta/sqrt(sum(w*delta^2))*sqrt(n*m)       #normalize dissimilarities
  
  #delta_plus <- ifelse(delta>=0,delta,0)  #delta decomposition (+)
  #delta_min <- ifelse(delta<=0,-delta,0)  #delta decomposition (-) (if all >0 --> complete 0)
  
  if (is.list(init)) {
    x <-init[[1]]                         #list as input structure
    y <-init[[2]]
  } else {
    e <- dhat^2
    e <- -0.5*(e-outer(rowSums(e)/m,colSums(e)/n,"+")+(sum(e)/(n*m)))
    #e <- e/sqrt(sum(e^2))*sqrt(n*m) 
    
    z <- svd(e,nu=p,nv=0)                 #SVD for e (pos. distances)
    x<-z$u                                #starting value for x
    y<-crossprod(e,x)                     #starting value for y
  }
  
  if (circle != "none"){
    r <- projCircle(x,y,x,y,circle=circle)
    x <- r$x
    y <- r$y
    wr <- rowSums(w)
    wc <- colSums(w)
    lambda <- 2*max(c(wr,wc))    
  }
  
  d <- distRect(x,y,reg)                  #n times m of reproduced diss
  coefOfVar <- function(x, w){            # Compute coefficient of variation
    av <- sum(x * w)/sum(w)
    va <- sum((x - av)^2 * w)/sum(w)
    return(va/av)
  }
  pstress <- function(dhat, d, w, omega, lambda, sumw){
    if (conditionality == "matrix"){
      g <- rep(0, 3)
      g[1] <- sum(w*(dhat - d)^2)/sumw
      av <- sum(dhat * w)/sumw
      va <- sum(dhat^2 * w) - sumw * av^2
      g[2] <- va + omega * sumw * av^2
      g[3] <- va^.5
    } else { ## conditionality == "row")
      pen <- rep(NA, n)
      for (i in 1:n){
        pen[i] <- 1/coefOfVar(dhat[i, ], w[i, ])
      }
    }
    return(list(pstress = g[1]^lambda * g[2]/g[3]^2, g = g, va = va, av = av))
    
  }
  
  ps   <- pstress(dhat, d, w, omega, lambda, sumw)   #pstress value
  lold <- ps$pstress

  
  
  #------------------- begin majorization -----------------------------
  repeat {
    if (circle == "none") {
      ww <- w 
      wr <- rowSums(ww)
      wc <- colSums(ww)
      
      v <- solve(diag(wc) + (1/m) - crossprod(ww, ww/wr)) - (1/m)
      
      b  <- w * dhat / d                 #B matrix
      br <- rowSums(b)                   #rows B
      bc <- colSums(b)                   #columns W
      
      xraw <- (br * x) - ( b %*% y)
      yraw <- (bc * y) - crossprod(b, x)
      
      y <- v %*% (yraw + crossprod(ww, xraw/wr))#x update 
      x <- (xraw + (ww %*% y))/wr                #y update
      
    } else {
      b  <- w*(1-dhat/d)                #B matrix
      br <- rowSums(b)                   #rows B
      bc <- colSums(b)                   #columns W
      xunc <- x - outer(br, rep(1/lambda, p), "*") * x + b %*% (y/lambda)
      yunc <- y - outer(bc, rep(1/lambda, p), "*") * y + t(b) %*% (x/lambda)
      r <- projCircle(xunc, yunc, x, y, circle = circle)
      x <- r$x
      y <- r$y
    }  
    
    d <- distRect(x,y,reg)             #compute distances (update)
    
    #lnew <- sum(w*(dhat - d)^2)/sumw     #compute stress
    
    # Update dhats
    
    dhat.old <- dhat
    if (conditionality == "matrix"){
      alpha2 <- 0.5 * lambda * ps$g[2]^.5 * ps$g[1]^(lambda/2 - 1)
      alpha3 <- 0.5 * ps$g[1]^(lambda/2) * ps$g[2]^-0.5
      g      <- ps$pstress^.5
      tau1   <- sum(w * dhat) / sumw
      tau2   <- min(1/(2 * min(dhat)), 1/eps)
      tau4   <- tau1 / ps$g[3]
      beta1  <- 1/sumw
      beta3  <- 1 + omega + 2 * omega * tau1 * tau2
      beta5  <- tau2 * tau4
      cm     <- alpha2 * beta1 + alpha3 * beta3 + g * beta5
      t1     <- tau2 * as.vector(dhat) - 0.5
      t2     <- as.vector(dhat)/(2 * ps$va^.5)
      b1     <- as.vector(d)/sumw
      b2     <- tau1 + omega * as.vector(dhat) + (2*omega*tau1) * t1
      b3     <- t2 + tau4 * t1
      ksi    <- (alpha2 * b1 + alpha3 * b2 + g * b3)/cm
      tt     <- transform(ksi, disobj[[1]], w = as.vector(w))
      dhat   <- matrix(tt$res, n, m)  ## dhat update
    } else { ## conditionality == "row"
      dhat <- matrix(0, n, m)
      for (i in 1:n) {
        dhat[i, ] <- transform(d[i, ], disobj[[i]], w = w[i, ])$res  ## dhat update
      }
    }
    
    ps   <- pstress(dhat, d, w, omega, lambda, sumw)   #pstress value
    lnew <- ps$pstress
    
    if (verbose) cat("Iteration: ", formatC(itel, digits=6, width=6), 
                     "   Stress:",  formatC(lnew, digits=6, width=12, format="f"),
                     "   Dif:",     formatC(lold - lnew, digits=6,width=12, format="f"),
                     "\n")
    
    if ( ( (lold-lnew) < eps & itel > 1) || (itel==itmax)) break() 
    
    lold <- lnew                       #update stress
    itel <- itel+1
  }
  #-------------------- end majorization --------------------------
  
  colnames(y) <- colnames(x) <- paste("D",1:(dim(y)[2]),sep="")
  rownames(x) <- rownames(diss) <- rownames(d) <- rnames
  
  # point stress 
  resmat <- as.matrix(d - diss)^2    #point stress
  spp.col <- colMeans(resmat, na.rm = TRUE)
  spp.col <- spp.col/sum(spp.col)*100
  spp.row <- rowMeans(resmat, na.rm = TRUE)
  spp.row <- spp.row/sum(spp.row)*100
  
  if (itel == itmax) warning("Iteration limit reached! Increase itmax argument!")
  
  ## stress normalization
  lnew <- sqrt(sum(w*(dhat - d)^2, na.rm = TRUE)/sumw)
  
  ## congruence coefficients
  diss0 <- diss
  diss0[is.na(diss0)] <- 0
  congnum <- diag(diss0 %*% t(d))
  congdenom <- sqrt(diag(diss0 %*% t(diss0)) * diag(d %*% t(d)))
  congvec <- congnum/congdenom
  
  #return configuration distances, row and column configurations, stress 
  result <- list(obsdiss = diss, confdist = d, dhat = dhat, iord = tt$iord.prim, conf.row = x, conf.col = y, stress = lnew, 
                 spp.row = spp.row, spp.col = spp.col, congvec = congvec, weightmat = w,
                 ndim = p, model = "Rectangular smacof", niter = itel, nind = n, nobj = m, call = match.call()) 
  class(result) <- c("smacofR", "prefscal")
  result 
}