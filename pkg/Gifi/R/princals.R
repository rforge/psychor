princals <- function (data, ndim = 2, ordinal = TRUE, ties = "s", knots = knotsQ(data), degrees = 2, copies = 1, 
                      missing = "m", active = TRUE, itmax = 1000, eps = 1e-6, seed = 123, verbose = FALSE)  {
    
    ## --- sanity checks
    names <- colnames(data, do.NULL = FALSE) 
    rnames <- rownames(data, do.NULL = FALSE)
    data_orig <- data
    data <- makeNumeric(data)            
    ## insert missing and ties match.arg
    ## --- end sanity checks
  
    aname <- deparse (substitute (data))
    nvars <- ncol(data)
    nobs <- nrow(data)
    g <- makeGifi(data = data, knots = knots, degrees = reshape (degrees, nvars), ordinal = reshape(ordinal, nvars),
                  sets =  1:nvars, copies = reshape (copies, nvars), ties = reshape (ties, nvars), missing = reshape (missing, nvars),
                  active = reshape (active, nvars), names = names)
    h <- gifiEngine(gifi = g, ndim = ndim, itmax = itmax, eps = eps, seed = seed, verbose = verbose)
    a <- v <- z <- d <- y <- o <- as.list (1:nvars)
    dsum <- matrix (0, ndim, ndim)
    for (j in 1:nvars) {
      jgifi <- h$xGifi[[j]][[1]]
      v[[j]] <- jgifi$transform
      a[[j]] <- jgifi$weights
      y[[j]] <- jgifi$scores
      z[[j]] <- jgifi$quantifications
      cy <- crossprod (y[[j]])
      dsum <- dsum + cy
      d[[j]] <- cy
      o[[j]] <- crossprod (h$x, v[[j]])
    }
    
    # return (structure(list(transform = v, rhat = corList(v), objectscores = h$x, scores = y, quantifications = z,
    #                        dmeasures = d, lambda = dsum/ncol (data), weights = a, loadings = o, ntel = h$ntel, f = h$f),
    #   class = "princals"
    # ))
    
   ## --- output cosmetics
   dnames <- paste0("D", 1:ndim)
   transform <- do.call(cbind, v); colnames(transform) <- names; rownames(transform) <- rnames
   rhat <- corList(v); rownames(rhat) <- colnames(rhat) <- names
   evals <- eigen(rhat)$values
   objectscores <- as.matrix(h$x); colnames(objectscores) <- dnames; rownames(objectscores) <- rnames
   scoremat <- sapply(y, function(xx) xx[,1]); colnames(scoremat) <- names; rownames(scoremat) <- rnames
   quantifications <- as.matrix(z); names(quantifications) <- names; quantifications <- lapply(quantifications, "colnames<-", dnames)
   dmeasures <- d; names(dmeasures) <- names; dmeasures <- lapply(dmeasures, "colnames<-", dnames); dmeasures <- lapply(dmeasures, "rownames<-", dnames)
   lambda <- dsum/ncol(data); rownames(lambda) <- colnames(lambda) <- dnames
   weights <- do.call(rbind, a); rownames(weights) <- names; colnames(weights) <- dnames 
   loadings <- do.call(cbind, o); rownames(loadings) <- dnames; colnames(loadings) <- names 
   ntel <- h$ntel
   f <- h$f
   
   res <- list(transform = transform, rhat = rhat, evals = evals, objectscores = objectscores, scoremat = scoremat, quantifications = quantifications,
               dmeasures = dmeasures, lambda = lambda, weights = weights, loadings = loadings, ntel = ntel, f = f, 
               data = data_orig, datanum = data, ndim = ndim, 
               call = match.call())
   class(res) <- "princals"
   return(res)
  }