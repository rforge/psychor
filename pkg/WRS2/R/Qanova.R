Qanova <- function(formula, data, q = 0.5, nboot = 600){
  #
  # Test global hypothesis that J independent groups
  # have equal medians.
  # Performs well when there are tied values.
  #
  # Basically, use pbadepth in conjunction with the Harrell--Davis
  # estimator.
  #
  #
  if (missing(data)) {
    mf <- model.frame(formula)
  } else {
    mf <- model.frame(formula, data)
  }
  cl <- match.call()
  
  x <- split(model.extract(mf, "response"), mf[,2])   
  
  op=3
  MC <- FALSE
  
  chkcar=NA
  for(j in 1:length(x))chkcar[j]=length(unique(x[[j]]))
  if(min(chkcar<20)) warning("Cardinality of sample space is less than 20 for one more groups. Type I error might not be controlled!")

  output <- pbadepth1(x, est = hd, q = q, allp = TRUE, SEED = FALSE, op = op, nboot = nboot, MC = MC, na.rm = TRUE)
  result <- list(psihat = output$psihat, p.value = output$p.value, contrasts = output$contrasts, call = cl)
  class(result) <- "qanova"
  result
}