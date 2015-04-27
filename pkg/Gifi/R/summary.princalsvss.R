summary.princalsvss <- function(object, ...) UseMethod("summary.princalsvss")

summary.princalsvss <- function(object, ...)
{
  vsstab <- data.frame(VSS1 = object$VSS1, VSS2 = object$VSS2, 
                  MAP = object$MAP, SQres = object$sqresid, Fit = object$fit)
  rownames(vsstab) <- paste0("Comp.", 1:length(object$VSS1))
  res <- list(vsstab = vsstab)
  class(res) <- "summary.princalsvss"
  res
}

