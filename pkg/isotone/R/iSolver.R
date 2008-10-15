# SILF Loss

# SILF Loss

iSolver<-function(x,a,extra) {
w<-extra$w; z<-extra$z; eps<-extra$eps; beta<-extra$beta
fobj<-function(x) {
    y<-abs(x-z)
    g<-((y-(1-beta)*eps)^2)/(4*beta*eps)
    g[which(y < (1-beta)*eps)]<-0
    ii<-which(y > (1+beta)*eps)
    g[ii]<-y[ii]-eps
    return(sum(w*g))
    }
gobj<-function(x) {
    y<-x-z
    g<-rep(0,length(y))
    g[which(y < -(1+beta)*eps)]<--1
    ii<-which((y > -(1+beta)*eps) & (y < -(1-beta)*eps))
    g[ii]<-(y[ii]+(1-beta)*eps)/(2*beta*eps)
    ii<-which((y > (1-beta)*eps) & (y < (1+beta)*eps))
    g[ii]<-(y[ii]-(1-beta)*eps)/(2*beta*eps)
    g[which(y > (1+beta)*eps)]<-1
    return(w*g)
    }
return(fSolver(x,a,list(fobj=fobj,gobj=gobj)))
}