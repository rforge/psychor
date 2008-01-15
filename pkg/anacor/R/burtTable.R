burtTable<-function(data, pf = FALSE) 
{
# burtTable makes a Burt matrix out of a data-frame or a
# profile frequency matrix

nvar<-dim(data)[2]
if(pf) nvar<-nvar-1
ff<-rep(1,dim(data)[1])
if (pf) ff<-data[,dim(data)[2]]
ncat<-sapply(1:nvar,function(j) length(table(data[,j])))
tcat<-sum(ncat)
first<-cumsum(c(1,ncat))[1:nvar]
burt<-matrix(0,tcat,tcat)
for (i in 1:nvar) {
	ii<-first[i]:(first[i]+ncat[i]-1); dd<-as.factor(data[,i])
	gi<-ifelse(outer(dd,levels(dd),"=="),1,0)
	if (i < nvar) {
		for (j in (i+1):nvar) {
			jj<-first[j]:(first[j]+ncat[j]-1)
                        dd<-as.factor(data[,j])
			gj<-ifelse(outer(dd,levels(dd),"=="),1,0)
			burt[jj,ii]<-t(burt[ii,jj]<-crossprod(gi,gj*ff))
			}
		}
	burt[ii,ii]<-crossprod(gi,gi*ff)
	}
return(burt)
}