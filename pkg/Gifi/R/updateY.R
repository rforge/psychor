`updateY` <-
function(dframe,x,y,active,rank,level,sets,verbose=0){
  nobj<-dim(x)[1]
  ndim<-dim(x)[2]
  nset<-length(sets)
  for (l in 1:nset) {                     ## number of sets (if no sets specified, sets = nvar)
  	indi <- sets[[l]]                     ## set index
    jndi <- indi[which(active[indi])]     ## index for active variables in set 
    if (length(jndi) == 0) next()         ## no active variables
  	
    ii <- which(!is.na(dframe[,jndi[1]])) ## non-NA values in active variables (index for objects)
    if (length(ii) == 0) next()           ## only NA
  	
    ss <- sumSet(dframe,nobj,ndim,y,jndi) ## expand catscores to data 0 where NA
    
    for (j in jndi) {                     ## runs over active variables in set
		  gg <- dframe[ii,j]                  ## original categories
      yy <- y[[j]]                        ## category scores for variable j
      d <- as.vector(table(gg))           ## category frequencies  
		  s1 <- sum((x[ii,]-ss[ii,])^2)       ## SSQ X - Gy
		  ss[ii,] <- ss[ii,]-yy[gg,]          ## just relevant for sets, otherwise 0
		  yc <- computeY(gg,x[ii,]-ss[ii,])   ## compute scores (just relevant for sets)
		  yy <- restrictY(d,yc,rank[j],level[j],verbose=verbose)$y  ## crucial!! optimal scaling
		  ss[ii,] <- ss[ii,]+yy[gg,]
		  s2 <- sum((x[ii,]-ss[ii,])^2)
		  y[[j]]<-yy		  	
   	}
  }
return(y)
}

