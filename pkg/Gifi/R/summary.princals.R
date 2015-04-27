summary.princals <- function(object, ...) UseMethod("summary.princals")

summary.princals <- function(object, ...)
{
  cumvar <- round(cumsum(object$vaf), 4)           ## Variance Accounted for
  evar <- round(object$vaf, 4)
  comptab <- rbind(evar, cumvar)
  colnames(comptab) <- paste("Comp.", 1:length(evar), sep = "")
  rownames(comptab) <- c("VAF (%)", "Cumulative VAF (%)")
  res <- list(vartab = comptab, loadings = object$loadings)
  class(res) <- "summary.princals"
  res
}

