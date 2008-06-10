`lineals` <-
function(data, itmax = 100, eps = 1e-6)
{
  
  m <- dim(data)[2]                                      #number of variables
  n <- dim(data)[1]                                      #number of observations
  r <- diag(m)
  t <- diag(m)
  fold <- Inf
  itel <- 1
  
  ncat <- sapply(1:m, function(j) length(table(data[,j])))  #number of categories for each variable
  ccat <- c(0,cumsum(ncat))
  
  y <- list()
  for (j in 1:m) y <- c(y,list(1:ncat[j]))              #list with category indices (inital score values)
  names(y) <- colnames(data)

  burt <- burtTable(data)                               #create Burt matrix (each variable is taken as categorical)
  nameslist <- apply(data, 2, unique)
  colnames(burt) <- rownames(burt) <- unlist(lapply(nameslist, sort))
  d <- diag(burt)                                       #diagonal of the Burt matrix

  for (j in 1:m) {                                      #compute category scores (normalized to y'burt y = 1)
    indj <- (ccat[j]+1):ccat[j+1]
    dj <- d[indj]
    y[[j]] <- y[[j]] - sum(dj*y[[j]])/n
    y[[j]] <- y[[j]]/sqrt(sum(dj*y[[j]]*y[[j]]))        #burt normalization
  }
  
  #--------- start lineals iterations (maximize f) ----------------
  repeat {
    f <- 0                                              #loss value
    for (j in 1:m) {
      indj <- (ccat[j]+1):ccat[j+1]
      yj <- y[[j]]
      for (l in 1:m) {
	indl <- (ccat[l]+1):ccat[l+1]
        dl <- d[indl]
        yl <- y[[l]]
	r[j,l] <- sum(burt[indj,indl]*outer(yj,yl))                  #correlation matrix
	c <- burt[indj,indl]%*%diag(1/pmax(1,dl))%*%burt[indl,indj]
	t[j,l] <- sum(c*outer(yj,yj))			             #correlation ratios
	f <- f+(t[j,l]-r[j,l]^2)                                     #loss update (cf. p. 448, de Leeuw 1988)
      }
    }

    for (j in 1:m) {                                                 #score updating  
	indj <- (ccat[j]+1):ccat[j+1]
        nc <- ncat[j]
        yj <- y[[j]]
	dj <- d[indj]
        c <- matrix(0,nc,nc)
	for (l in 1:m) {
	  if (j != l) {
	    indl <- (ccat[l]+1):ccat[l+1]
            dl <- d[indl]
            yl <- y[[l]]
	    u <- burt[indj,indl]%*%(diag(1/pmax(1,dl)) - 2*outer(yl,yl))%*%burt[indl,indj] 
	    c <- c+u 
	  }
	}
	e <- eigen(c/sqrt(outer(dj,dj)))
	y[[j]] <- e$vectors[,nc]/sqrt(dj)                              #scores update
        #FIXME!!! incorporate optimal scaling
     }
     if (((fold-f) < eps) || (itel == itmax)) break
     itel <- itel+1
     fold <- f
  }

  if (itel == itmax) warning("Maximum iteration limit reached!")

  dummy.mat <- as.matrix(expandFrame(data))
  scorevec <- unlist(y)
  scoremat <- t(apply(dummy.mat, 1, function(xx) scorevec[which(xx == 1)]))
  colnames(scoremat) <- colnames(data)
  rownames(r) <- colnames(r) <- colnames(t) <- rownames(t) <- colnames(scoremat)
  for (i in 1:length(y)) names(y[[i]]) <- unique(data[,i])
  
  result <- list(loss = f, catscores = y, cormat = r, cor.rat = t, indmat = dummy.mat, 
  scoremat = scoremat, data = data, burtmat = burt, niter = itel, call = match.call())
  class(result) <- "aspect"
  result
}

