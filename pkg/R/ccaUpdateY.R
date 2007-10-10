`ccaUpdateY` <-
function(dframe,x,y,totsum,active,rank,level,sets){
  nobj<-dim(x)[1]; ndim<-dim(x)[2]; nset<-length(sets)
  for (l in 1:nset) {
  	indi<-sets[[l]]
  	ss<-sumSet(dframe,nobj,ndim,y,indi,active)
  	for (j in indi) {
  		if (active[j]){
  			gg<-dframe[,j]; yy<-y[[j]]; ii<-which(!is.na(gg)); d<-as.vector(table(gg))
  			ss[ii,]<-ss[ii,]-yy[gg[ii],]
  			ycen<-computeY(gg[ii],x[ii,]-ss[ii,])
  			yy<-restrictY(d,ycen,rank[j],level[j])$y
  			ss[ii,]<-ss[ii,]+yy[gg[ii],]
  			y[[j]]<-yy
  			}
  		}
  	totsum<-totsum+ss
  	}
  return(list(y=y,totsum=totsum))
}

