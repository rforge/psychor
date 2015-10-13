function (delta, constraint = "linear", external, ndim = 2, type = c("ratio", "interval", "ordinal", "mspline"), 
          weightmat = NULL, init = NULL, 
          ties = "primary", verbose = FALSE, modulus = 1, itmax = 1000, 
          eps = 1e-06, spline.intKnots = 4, spline.degree = 2, 
          constraint.type = c("ratio", "interval", "ordinal", "spline", "mspline"), constraint.ties = "primary", 
          constraint.spline.intKnots = 2, constraint.spline.degree = 2) 
{
  type <- match.arg(type, c("ratio", "interval", "ordinal", "mspline"), several.ok = FALSE)
  ties <- match.arg(ties, c("primary", "secondary", "tertiary"), several.ok = FALSE)
  constraint.type <- match.arg(constraint.type, c("ratio", "interval", "ordinal", "spline", "mspline"), several.ok = FALSE)
  constraint.ties <- match.arg(constraint.ties, c("primary", "secondary", "tertiary"), several.ok = FALSE)
  diss <- delta
  if ((is.matrix(diss)) || (is.data.frame(diss))) {
    diss <- strucprep(diss)
    attr(diss, "Labels") <- rownames(delta)
  }
  checkdiss(diss)
  p <- ndim
  n <- attr(diss, "Size")
  if (p > (n - 1)) 
    stop("Maximum number of dimensions is n-1!")
  nn <- n * (n - 1)/2
  m <- length(diss)
  startconf <- init
  if (!is.null(startconf)) 
    startconf <- as.matrix(init)
  xstart <- startconf
  if (is.null(attr(diss, "Labels"))) 
    attr(diss, "Labels") <- paste(1:n)
  if (is.data.frame(external)) 
    external <- as.matrix(external)
  if (!is.list(external)) {
    if (ncol(external) < p) 
      stop("Number of external variables can not be smaller than the number of MDS dimensions!")
  }
  if (is.matrix(external)) {
    extvars <- list()
    for (s in 1:ncol(external)) {
      constraint.trans <- constraint.type
      if (constraint.trans == "ratio") {
        constraint.trans <- "none"
      }
      else if (constraint.trans == "ordinal" & ties == "primary") {
        constraint.trans <- "ordinalp"
      }
      else if (constraint.trans == "ordinal" & ties == "secondary") {
        constraint.trans <- "ordinals"
      }
      else if (constraint.trans == "ordinal" & ties == "tertiary") {
        constraint.trans <- "ordinalt"
      }
      else if (constraint.trans == "spline") {
        constraint.trans <- "spline"
      }
      else if (constraint.trans == "mspline") {
        constraint.trans <- "mspline"
      }
      extvars[[s]] <- transPrep(external[, s] - mean(external[, s], na.rm = TRUE), trans = constraint.trans, 
                                spline.intKnots = constraint.spline.intKnots, 
                                spline.degree = constraint.spline.degree, missing = "multiple")
      external[, s] <- extvars[[s]]$xInit - mean(extvars[[s]]$xInit)
    }
  }
  simpcirc <- FALSE
  if (is.list(external)) {
    if (external[[1]] == "simplex") {
      d2 <- external[[2]]
      if (d2 >= n) 
        stop("Simplex dimension must be < n!")
      external <- diag(1, n)[, 1:d2]
      external[lower.tri(external)] <- 1
    }
    if (external[[1]] == "circumplex") {
      d2 <- external[[2]]
      if (d2 >= n) 
        stop("Circumplex dimension must be <= n!")
      k1 <- external[[3]]
      k2 <- external[[4]]
      if (k2 <= k1) 
        stop("k2 must be > k1")
      external <- matrix(0, nrow = n, ncol = d2)
      ind.all <- expand.grid(1:n, 1:d2)
      inddiff <- apply(ind.all, 1, function(xx) abs(diff(xx)))
      ind.good <- which((inddiff >= k1) + (inddiff <= k2) == 2)
      el1 <- as.matrix(ind.all[ind.good, ])
      external[el1] <- 1
    }
    simpcirc <- TRUE
  }
  K <- dim(external)[2]
  if (is.null(weightmat)) {
    wgths <- initWeights(diss)
  }
  else wgths <- as.dist(weightmat)
  trans <- type
  if (trans == "ratio") {
    trans <- "none"
  }
  else if (trans == "ordinal" & ties == "primary") {
    trans <- "ordinalp"
  }
  else if (trans == "ordinal" & ties == "secondary") {
    trans <- "ordinals"
  }
  else if (trans == "ordinal" & ties == "tertiary") {
    trans <- "ordinalt"
  }
  else if (trans == "spline") {
    trans <- "mspline"
  }
  disobj <- transPrep(diss, trans = trans, spline.intKnots = spline.intKnots, spline.degree = spline.degree)
  dhat <- normDissN(diss, wgths, 1)
  dhat[is.na(dhat)] <- 1
  w <- vmat(wgths)
  v <- myGenInv(w)
  itel <- 1
  if (!is.function(constraint)) {
    if (constraint == "linear") {
      constrfun <- function(x, w, external) {
        external %*% solve(crossprod(external, w %*% external), crossprod(external, w %*% x))
      }
      if (is.null(xstart)) 
        xstart <- matrix(rnorm(n * p), n, p)
    }
    if (constraint == "diagonal") {
      constrfun <- function(x, w, external) {
        return(external %*% diag(colSums(external * (w %*% x))/colSums(external * (w %*% external))))
      }
      if (is.null(xstart)) 
        xstart <- matrix(rnorm(n * K), n, K)
    }
    if (constraint == "unique") {
      constrfun <- function(x, w, external) {
        n <- dim(x)[1]
        p <- dim(x)[2] - n
        return(cbind(x[, 1:p], diag(diag(w %*% x[, p + (1:n)])/diag(w))))
      }
      if (is.null(xstart)) 
        xstart <- cbind(matrix(rnorm(n * p), n, p), diag(1, n))
    }
  }
  else {
    constrfun <- constraint
    if (is.null(xstart)) 
      stop("Starting configuration must be specified!")
  }
  if (constraint %in% c("linear", "diagonal") & !simpcirc) {
    ncol.ext <- ncol(external)
    if (constraint == "linear") {
      C <- svd(diag(ncol.ext) - 1/ncol.ext)$u[, 1:ndim]
    }
    else if (constraint == "diagonal") {
      C <- diag(ncol.ext)
    }
    x.unc <- xstart
    for (s in 1:ncol.ext) {
      target <- x.unc %*% C[s, ]
      external[, s] <- target
      tt.plus <- transform(target, extvars[[s]], normq = 0)
      tt.min <- transform(-target, extvars[[s]], normq = 0)
      if (sum((tt.plus$res - target)^2) < sum(((tt.min$res + target))^2)) {
        external[, s] <- tt.plus$res
      }
      else {
        external[, s] <- tt.min$res
      }
      x.unc <- x.unc - outer(external[, s], C[s, ])
    }
    updext.result <- updext(xstart, w, external, extvars, constraint)
    external <- updext.result$external
  }
  x <- constrfun(xstart, w, external)
  d <- dist(x)
  lb <- sum(wgths * d * dhat)/sum(wgths * d^2)
  x <- lb * x
  d <- lb * d
  sold <- sum(wgths * (dhat - d)^2)/nn
  repeat {
    b <- bmat(dhat, wgths, d)
    y <- v %*% b %*% x
    if (constraint %in% c("linear", "diagonal") & !simpcirc) {
      updext.result <- updext(x, w, external, extvars, constraint)
      external <- updext.result$external
    }
    y <- constrfun(y, w, external)
    e <- dist(y)
    ssma <- sum(wgths * (dhat - e)^2)
    dhat2 <- transform(e, disobj, w = wgths, normq = nn)
    dhat <- dhat2$res
    snon <- sum(wgths * (dhat - e)^2)/nn
    if (verbose) 
      cat("Iteration: ", formatC(itel, width = 3, format = "d"), 
          " Stress (raw): ", formatC(c(snon), digits = 8, width = 10, format = "f"), " Difference: ", 
          formatC(c(sold - snon), digits = 8, width = 10, 
                  format = "f"), "\n")
    if (((sold - snon) < eps) || (itel == itmax)) 
      (break)()
    x <- y
    d <- e
    sold <- snon
    itel <- itel + 1
  }
  stress <- sqrt(snon)
  if (any(is.na(y))) {
    csy <- colSums(y)
    ind <- which(is.na(csy))
    y <- y[, -ind]
  }
  colnames(y) <- paste("D", 1:(dim(y)[2]), sep = "")
  rownames(y) <- labels(diss)
  dhat <- structure(dhat, Size = n, call = quote(as.dist.default(m = b)), class = "dist", Diag = FALSE, Upper = FALSE)
  attr(dhat, "Labels") <- labels(diss)
  attr(e, "Labels") <- labels(diss)
  dhat[is.na(diss)] <- NA
  confdiss <- normDissN(e, wgths, 1)
  spoint <- spp(dhat, confdiss, wgths)
  rss <- sum(spoint$resmat[lower.tri(spoint$resmat)])
  if ((constraint == "diagonal") && (!simpcirc)) {
    if (p != ncol(y)) {
      warning("Number of dimensions is set equal to number of external variables!")
      p <- ncol(y)
    }
  }
  if (simpcirc) {
    if (p != ncol(y)) {
      warning("Number of dimensions is set equal to dimension defined by simplex/circumplex")
      p <- ncol(y)
    }
  }
  if (itel == itmax) 
    warning("Iteration limit reached! Increase itmax argument!")
  Z <- as.matrix(external)
  X <- y
  C <- solve(t(Z) %*% Z) %*% t(Z) %*% X
  if (constraint %in% c("linear", "diagonal") && !simpcirc) {
    for (s in 1:ncol(external)) {
      extvars[[s]]$iord.prim <- updext.result$iord.prim[[s]]
      extvars[[s]]$final <- external[, s]
      extvars[[s]]$c <- C[s, ]
    }
  }
  else {
    extvars = NULL
  }
  result <- list(delta = diss, dhat = dhat, confdiss = confdiss, 
                 conf = y, C = C, stress = stress, spp = spoint$spp, ndim = p, 
                 iord = dhat2$iord.prim, extvars = extvars, external = external, 
                 weightmat = wgths, resmat = spoint$resmat, rss = rss, 
                 init = xstart, model = "SMACOF constraint", niter = itel, 
                 nobj = n, type = type, call = match.call())
  class(result) <- c("smacofB", "smacof")
  result
}