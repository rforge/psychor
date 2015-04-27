## monotone spline transformation
monspline <- function (y, X) {
  
  a <- difmat(nrow(X))
  pa <- - lm.fit(X, t(a))$fitted.values
  pb <- drop(lm.fit(X, y)$fitted.values)
  lb <- nnls(pa, y)$x
  Xb <- pb - drop(pa %*% lb)
  return(list(yhat = Xb, RSS = sum((y - Xb) ^ 2)))
}