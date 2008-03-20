#smacof with linear constraints on the configuration (de Leeuw & Heiser, 1980; Borg & Groenen, p. 236)

smacofConstraint <- function(delta, constraint = "linear", external, ndim = 2, weightmat = NULL, metric = TRUE,
                             ties = "primary", verbose = FALSE, modulus = 1, itmax = 1000, eps = 1e-6)
{
# diss ... dissimilarity matrix
# constraint ... either "linear", "unique", "diagonal", or a user-specified function
# external ... external data for X-decomposition (Z in paper)
# wghts ... weight structure. if not specified, weights is 1-structure
# ndim ... number of dimensions
# metric ... if TRUE, metric MDS, if FALSE, non-metric
# ties ... ties for pava (primary, secondary, tertiary)
# modulus ... modulus for nonmetric update
# itmax ... maximum number of iterations
# eps ... change in loss function

  diss <- delta
  if ((is.matrix(diss)) || (is.data.frame(diss))) {
    diss <- strucprep(diss)  #if data are provided as dissimilarity matrix
    attr(diss, "Labels") <- rownames(delta)
  }
  p <- ndim                                     
  n <- attr(diss,"Size")
  m <- length(diss)
  
  if (is.null(attr(diss, "Labels"))) attr(diss, "Labels") <- paste(1:n)

  external <- as.matrix(external)
  
  if (is.null(weightmat)) {
    wgths <- initWeights(diss)
  }  else  wgths <- weightmat
  
 
  dhat <- normDiss(diss,wgths)                  #normalize dissimilarities
  w <- vmat(wgths)                              #matrix V
  v <- myGenInv(w)                              #Moore-Penrose inverse
  itel <- 1

  #----------- pre-specified functions for constraints -----------
  # linear constraint (de Leeuw & Heiser, 1980, p.515), X=ZC 
 if (!is.function(constraint)) { 
  if (constraint == "linear") {
    constrfun <-function(x,w,external) {
      return(external%*%solve(crossprod(external,w%*%external),crossprod(external,w%*%x)))
    }
  }
  
  # C restricted to be diagonal   
  if (constraint == "diagonal") {
    constrfun <- function(x,w,external) {
      return(external%*%diag(colSums(external*(w%*%x))/colSums(external*(w%*%external))))  
    }
  }
    
  # X with uniqueness coordinates  
  if (constraint == "unique") {
    constrfun <- function(x,w,external) {
      n <- dim(x)[1]
      p <- dim(x)[2]-n 
      return(cbind(x[,1:p],diag(diag(w%*%x[,p+(1:n)])/diag(w))))
    }
  }
 } else {   # user-specified
  constrfun <- constraint
 }
    
  #---------- end pre-specified functions for constraints -------
  
  x <- constrfun(matrix(rnorm(n*p),n,p),w,external)    #compute X   
  d <- dist(x)                                         #distances X
  lb <- sum(wgths*d*dhat)/sum(wgths*d^2)               #denominator: normalization tr(X'VX) 
  x <- lb*x                                            #modify x with lb-factor
  d <- lb*d                                            #modify d with lb-factor
  sold <- sum(wgths*(dhat-d)^2)                        #initial stress

  #------- begin majorization -----------
  repeat {                                             #majorization iterations
	b <- bmat(dhat,wgths,d)                        #B matrix
        y <-v%*%b%*%x                                  #Y computation
	y <- constrfun(y,w,external)                   #update Y with corresponding constraints            
        e <- dist(y)                                   #Y distances
	ssma <- sum(wgths*(dhat-e)^2)                  #new stress

	if (!metric) {                                 #nonmetric MDS
	    if ((itel%%modulus) == 0) {
			if (ties=="primary") daux <- monregP(diss,e,wgths)
			if (ties=="secondary") daux <- monregS(diss,e,wgths)
			if (ties=="tertiary") daux <- monregT(diss,e,wgths)
			dhat<-normDiss(daux,wgths)
			}
	}
	snon <- sum(wgths*(dhat-e)^2)                  #nonmetric stress
        
	if (verbose) cat("Iteration: ",formatC(itel,width=3, format="d")," Stress: ",
		formatC(c(sold,ssma,snon),digits=8,width=12,format="f"),"\n")

	if (((sold-snon)<eps) || (itel == itmax)) break()   #convergence 

        x <- y                                         #updates
        d <- e
        sold <- snon
        itel <- itel+1	
  }
  #------- end majorization -----------
 
  if (metric) snon <- NULL          #no non-metric stress
  if (!metric) ssma <- NULL
 
  colnames(y) <- paste("D",1:(dim(y)[2]),sep="")
  rownames(y) <- labels(diss)
  dhat <- structure(dhat, Size = n, call = quote(as.dist.default(m=b)), class = "dist", Diag = FALSE, Upper = FALSE) 
  attr(dhat, "Labels") <- labels(diss)
  attr(e, "Labels") <- labels(diss)

#return configurations, configuration distances, normalized observed distances 
  result <- list(obsdiss = dhat, confdiss = e, conf = y, stress.m = ssma, stress.nm = snon,
               ndim = p, model = "SMACOF constraint", niter = itel, nobj = n) 
  class(result) <- c("smacofB","smacof")
  result 
}

