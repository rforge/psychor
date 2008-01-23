#SMACOF for individual differences (list of dissimilarity matrices)

smacofIndDiff <- function(delta, weightmat = NULL, ndim = 2, init = NULL, metric = TRUE,
                          ties = "primary", constraint = NULL, verbose = FALSE, modulus = 1,
                          itmax = 100, eps = 1e-6)
  
# delta ... list of input objects: either of class dist() or a symmetric matrix
# contstraint ... either NULL, "ident", "diag", "idio"
  

{
  diss <- delta
  p <- ndim
  wgths <- weightmat
  constr <- constraint

  #insert lapply for transforming dist structure into matrix
  
  if (!is.list(diss)) diss <- list(diss)
  if (is.null(weightmat)) wghts <- initWeights(diss)
  if (!is.list(wgths)) wgths <- list(wgths)

  
  n <- attr(diss[[1]],"Size")
  m <- length(diss)
  itel <- 1
      
  dr <- list()
  wr <- list()
  vr <- list()
  dh <- list()
      
  for (j in 1:m) {
	wr<-appendList(wr,vmat(wgths[[j]]))
	vr<-appendList(vr,myGenInv(wr[[j]]))
	dh<-appendList(dh,normDiss(diss[[j]],wgths[[j]]))
	}
  xr <-list()
  sold <- sf1 <- sf2 <- 0

  if (is.null(init)) {
    aconf<-torgerson(sumList(diss),p)
  } else aconf<-matrix(rnorm(n*p),n,p)

  
bconf<-repList(diag(p),m)
for (j in 1:m) {
	xr[[j]]<-aconf%*%bconf[[j]]
	dr[[j]]<-dist(xr[[j]]); sf1<-sf1+sum(wgths[[j]]*dr[[j]]*dh[[j]])
	sf2<-sf2+sum(wgths[[j]]*dr[[j]]^2)
	}
lb<-sf1/sf2
for (j in 1:m) {
	aconf<-lb*aconf
	xr[[j]]<-lb*xr[[j]]; dr[[j]]<-lb*dr[[j]]
	sold<-sold+sum(wgths[[j]]*(dh[[j]]-dr[[j]])^2)
	}
repeat {
    br<-list(); yr<-list(); er<-list(); sunc<-0
	for (j in 1:m) {
		br<-appendList(br,bmat(dh[[j]],wgths[[j]],dr[[j]]))
		yr<-appendList(yr,vr[[j]]%*%br[[j]]%*%xr[[j]])
		er<-appendList(er,dist(yr[[j]]))
		sunc<-sunc+sum(wgths[[j]]*(dh[[j]]-er[[j]])^2)
		}
	scon<-sunc

        if (!is.null(constr))
        {
		scon<-0; er<-list()
		if (constr=="ident") {
			z<-matrix(0,n,p); u<-matrix(0,n,n)
			for (j in 1:m) {
				z<-z+wr[[j]]%*%yr[[j]]
				u<-u+wr[[j]]
				}
			aconf<-myGenInv(u)%*%z; yr<-repList(aconf,m)
			}
		if (constr=="diag") {
			aux0<-matrix(0,n,p)
			for (j in 1:m) {
				aux1<-diag(crossprod(aconf,wr[[j]]%*%yr[[j]]))
				aux2<-diag(crossprod(aconf,wr[[j]]%*%aconf))
				bconf[[j]]<-diag(aux1/aux2)
				aux0<-aux0+(wr[[j]]%*%yr[[j]]%*%bconf[[j]])		    			    
				}
			for (s in 1:p) {
				aux1<-matrix(0,n,n)
				for (j in 1:m) 
					aux1<-aux1+(bconf[[j]][s,s]^2)*wr[[j]]
				aconf[,s]<-myGenInv(aux1)%*%aux0[,s]
				}
			for (j in 1:m)
				yr[[j]]<-aconf%*%bconf[[j]]
			}
		if (constr=="idio") {
			aux0<-matrix(0,n,p); auxk<-matrix(0,(n*p),(n*p))
			for (j in 1:m) {
				aux1<-crossprod(aconf,wr[[j]]%*%yr[[j]])
				aux2<-crossprod(aconf,wr[[j]]%*%aconf)
				auxb<-solve(aux2,aux1); bconf[[j]]<-auxb
				auxc<-crossprod(t(auxb))
				aux0<-aux0+(wr[[j]]%*%yr[[j]]%*%t(auxb))	
				auxk<-auxk+kronecker(auxc,wr[[j]])    			    
				}
			auxv<-kronecker(diag(p),matrix((1/n),n,n))
		    aconf<-matrix(solve(auxk+auxv,as.vector(aux0)),n,p)
			for (j in 1:m)
				yr[[j]]<-aconf%*%bconf[[j]]
			}
		for (j in 1:m) {
			er<-appendList(er,dist(yr[[j]]))
			scon<-scon+sum(wgths[[j]]*(dh[[j]]-er[[j]])^2)
			}
		}
	snon<-scon
	if (!metric) {
	    if ((itel%%modulus) == 0) {
			snon<-0; dh<-list()
			for(j in 1:m) {
				ds<-diss[[j]]; es<-er[[j]]; ws<-wgths[[j]]
				if (ties=="primary") do<-monregP(ds,es,ws)
				if (ties=="secondary") do<-monregS(ds,es,ws)
				if (ties=="tertiary") do<-monregT(ds,es,ws)
				dh<-appendList(dh,normDiss(do,ws))
				snon<-snon+sum(ws*(dh[[j]]-es)^2)
				}
			}
		}
	if (verbose) cat("Iteration: ",formatC(itel,width=3, format="d")," Stress: ",
		formatC(c(sold,sunc,scon,snon),digits=8,width=12,format="f"),"\n")
	if (((sold-snon)<eps) || (itel == itmax)) break()
	xr<-yr; dr<-er; sold<-snon; itel<-itel+1
	}
list(dh=dh,dr=er,xr=yr,a=aconf,b=bconf)
}

