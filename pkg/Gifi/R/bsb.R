## creates a bspline basis
bsb <- function (x, knots, order = 3) {
  
  if (length(knots) == 1) knots <- quantile(x, probs = seq(0, 1, 1/knots))
  if (min(x) < min(knots)) stop("B-spline argument smaller than smallest knot.")
  if (max(x) > max(knots)) stop("B-spline argument larger than largest knot.")
  
  nrw <- length(x)
  nbr <- length(knots)
  nrs <- order + nbr - 2
  g <- matrix (0, nrw, nrs)
  for (i in 1:nrw)
    g[i,] <- .C("BSPLINE",as.double(x[i]),as.integer(order),as.integer(nbr),as.double(knots),results = as.double(rep(0.0,nrs)))$results
  return(list(X = g[,which(colSums(g) > 0)], knots = knots))
}
