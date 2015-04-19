summary.princals <- function(object, ...) UseMethod("summary.princals")

summary.princals <- function(object, ...)
{
  cumvar <- round(cumsum(object$expvar), 4)
  evar <- round(object$expvar, 4)
  comptab <- rbind(evar, cumvar)
  colnames(comptab) <- paste("Comp.", 1:length(evar), sep = "")
  rownames(comptab) <- c("Explained Variance (%)", "Cumulative Variance (%)")
  res <- list(vartab = comptab)
  class(res) <- "summary.princals"
  res
}

