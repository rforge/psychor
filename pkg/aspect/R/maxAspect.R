`maxAspect` <-
function(data, aspect = "aspectSum", itmax = 100, eps=1e-6, extra = 1) 
{
# aspect ... function names: either "aspectSum", "aspectAbs", "aspectSMC", "aspectSumSMC", "aspectEigen", "aspectDeterminant" or a user-specified function

  m <- dim(data)[2]
  n <- dim(data)[1]
  r <- diag(m)
  fold <- -Inf
  itel<-1

  ncat <- sapply(1:m,function(j) length(table(data[,j])))
  ccat <- c(0,cumsum(ncat))
  y <- list()
  for (j in 1:m) y<-c(y,list(1:ncat[j]))
  names(y) <- colnames(data)

  burt <- burtTable(data)                                      #compute Burt matrix
  nameslist <- apply(data, 2, unique)
  colnames(burt) <- rownames(burt) <- unlist(lapply(nameslist, sort))
  d <- diag(burt)

  for (j in 1:m) {                                             #initial scaling of y such that y' burt y = 1
    indj <-(ccat[j]+1):ccat[j+1]
    dj <- d[indj]
    y[[j]] <- y[[j]]-sum(dj*y[[j]])/n
    y[[j]] <- y[[j]]/sqrt(sum(dj*y[[j]]*y[[j]])/n)
  }

  #----------------- aspect string/function check --------------
  if (!is.function(aspect)) {
    if (aspect == "aspectSum") aspectfun <- aspectSum
    if (aspect == "aspectAbs") aspectfun <- aspectAbs
    if (aspect == "aspectSMC") aspectfun <- aspectSMC
    if (aspect == "aspectSumSMC") aspectfun <- aspectSumSMC
    if (aspect == "aspectEigen") aspectfun <- aspectEigen
    if (aspect == "aspectDeterminant") aspectfun <- aspectDeterminant
  } else {
    aspectfun <- aspect
  }

  #---------------------- end aspect check ---------------------
  
  #-------------------------------------------------------------
  repeat {
    for (j in 1:m) {
       indj <- (ccat[j]+1):ccat[j+1]
	for (l in 1:m) {
	  indl <- (ccat[l]+1):ccat[l+1]
	  r[j,l] <- sum(y[[j]]*(burt[indj,indl]%*%y[[l]]))/n
	}
    }

    #FIXME!!! implement aspect as function
    a <- aspectfun(r, extra)                           #call aspect as a function of the correlation matrix r (and extra)
    f <- a$f                                        #value of the aspect function
    g <- a$g

    #print(f)

    for (j in 1:m) {
      indj<-(ccat[j]+1):ccat[j+1]; y[[j]]<-rep(0,ncat[j]); dj<-d[indj]
      for (l in 1:m) {
	indl<-(ccat[l]+1):ccat[l+1]
	if (j != l) y[[j]] <- y[[j]]+g[j,l]*burt[indj,indl]%*%y[[l]]
      }
      y[[j]] <- y[[j]]/dj
      y[[j]] <- y[[j]]-sum(dj*y[[j]])/n
      y[[j]] <- y[[j]]/sqrt(sum(dj*y[[j]]*y[[j]])/n)
    }
    if (((f-fold) < eps) || (itel == itmax)) break
    itel<-itel+1
    fold<-f
  }

  #returns: function value (total discrepancy), scores, correlation matrix of the scores, eigenvalues (of r)
  if (itel == itmax) warning("Maximum iteration limit reached!")

  dummy.mat <- as.matrix(expandFrame(data))
  scorevec <- unlist(y)
  scoremat <- t(apply(dummy.mat, 1, function(xx) scorevec[which(xx == 1)]))
  colnames(scoremat) <- colnames(data)
  rownames(r) <- colnames(r) <- colnames(scoremat)
  for (i in 1:length(y)) {
    rownames(y[[i]]) <- unique(data[,i])
    colnames(y[[i]]) <- "score"
  }
 
  
  result <- list(loss = f, catscores = y, cormat = r, eigencor = eigen(r,only.values=TRUE)$values, indmat = dummy.mat, scoremat = scoremat, data = data, burtmat = burt, niter = itel, call = match.call())
  class(result) <- "aspect"
  result
}

