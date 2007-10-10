`pcaUpdateY` <-
function(dframe,x,y,totsum,active,rank,level,nset){
  for (j in 1:nset){
  	if (active[j]){
  		gg<-dframe[,j]; ycen<-computeY(gg,x); d<-as.vector(table(gg)); ii<-which(!is.na(gg))
  		yy<-restrictY(d,ycen,rank[j],level[j])$y
  		y[[j]]<-yy
  		totsum[ii,]<-totsum[ii,]+yy[gg[ii],]
  		}
  	}
  return(list(y=y,totsum=totsum))
}

