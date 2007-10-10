`homals` <-
function(dframe,   # data (in data-frame)
sets=0,              	# list of vectors of indices
	ndim=2,              	# dimensionality (default 2)
	active=TRUE,            # which variables are active (single TRUE means all)
  rank=ndim,           	# which quantification ranks (default all ndim)
  level = "nominal",
  eps1=-Inf,           	# iteration precision eigenvalues (default 1e-6)
	eps2=1e-6,           	# iteration precision eigenvectors (default 1e-6)
	itermax=100         	# maximum number of iterations (default 100)
)

{

name <- deparse(substitute(dframe))
nobj <- dim(dframe)[1]                         #number of objects
nvar <- dim(dframe)[2]                         #number of variables
if ((length(level) > 1) && (length(level) != nvar)) 
  stop ("Level vector does not match number of variables!")


iter <- 0
pops <- 0

for (j in 1:nvar) {
	dframe[, j] <- as.factor(dframe[, j])
    levels(dframe[, j])<-sort(levels(dframe[, j]))
	}

if (length(sets) == 1) sets <- lapply(1:nvar,"c")              #variable sets
pca <- (max(sapply(sets,length)) <= 1)

nset <- length(sets)                                           #number of sets
if (pca) { mis <- apply(dframe,1,function (x) sum(ifelse(is.na(x),0,1)))
} else {
   mis <- apply(sapply(1:nset,
		function(s) apply(cbind(dframe[,sets[[s]]]),1,
		function (x) prod(ifelse(is.na(x),0,1)))),1,sum)
}

vname <- attr(dframe,"names")
rname <- attr(dframe,"row.names")

if (length(active)==1) active<-rep(active,nvar)

if (length(rank)==1) rank <- rep(rank,nvar)                  #quantification rank
if (length(level)==1) level<-rep(level,nvar)                 #nominal (quantification) level

for (j in 1:nvar) {
  k<-length(levels(dframe[,j]))
	if (rank[j] > min(ndim,k-1)) rank[j]<-min(ndim,k-1)
	}
	
x <- cbind(orthogonalPolynomials(mis,1:nobj,ndim))
x <- normX(centerX(x,mis),mis)$q
y <- lapply(1:nvar, function(j) computeY(dframe[,j],x))

#-----------------------------main computation--------------------------------
repeat {
	iter <- iter+1
	totsum <- array(0.0,dim(x))
	if (pca) {
    up.y <- pcaUpdateY(dframe,x,y,totsum,active,rank,level,nset)
  } else {
    up.y <- ccaUpdateY(dframe,x,y,totsum,active,rank,level,sets)
  }
	y <- up.y$y
  totsum <- up.y$totsum
	qv <- normX(centerX((1/mis)*totsum,mis),mis)
	z <- qv$q
  r<-qv$r
  ops <- sum(r)
  aps <- sum(La.svd(crossprod(x,mis*z),0,0)$d)/ndim
  if (iter == itermax) {
		warning("maximum number of iterations reached")
		break()
		}
  if (pops > ops) {
		warning(cat("loss function increases in iteration ",iter,"\n"))
		break()
		}
  if (((ops - pops) < eps1) || ((1.0 - aps) < eps2)) break()
		else {x<-z; pops<-ops}
}

ylist<-alist<-clist<-ulist<-NULL
for (j in 1:nvar) {
  gg<-dframe[,j]; c<-computeY(gg,z); d<-as.vector(table(gg))
  lst<-restrictY(d,c,rank[j],level[j])
  y<-lst$y; a<-lst$a; u<-lst$z
  ylist<-c(ylist,list(y)); alist<-c(alist,list(a)); clist<-c(clist,list(c)); ulist<-c(ulist,list(u))
}

#--------------------------preparing/labeling output----------------------------
dimlab <- paste("D", 1:ndim, sep = "")
for (i in 1:nvar) {
  rownames(ylist[[i]]) <- rownames(ulist[[i]]) <- rownames(clist[[i]])
  rownames(alist[[i]]) <- paste(1:dim(alist[[i]])[1])
  colnames(clist[[i]]) <- colnames(ylist[[i]]) <- colnames(alist[[i]]) <- dimlab
  colnames(ulist[[i]]) <- paste(1:dim(ulist[[i]])[2])
}
names(ylist) <- names(ulist) <- names(clist) <- names(alist) <- colnames(dframe)
rownames(z) <- rownames(dframe)
colnames(z) <- dimlab
#alist.t <- lapply(alist,t)
#--------------------------end preparing/labeling output------------------------

result <- list(datname = name, dframe = dframe, ndim = ndim, niter = iter, level = level, 
               eigenvalues = r, gain = c(ops, aps), rank.vec = rank,
               scores = z, rank.cat = ylist, cat.centroids = clist,
               cat.loadings = alist, low.rank = ulist)
class(result) <- "homals"
result
}

