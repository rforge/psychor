lfSolver<-function(x,a,extra) {
    w<-extra$w; z<-extra$z; n<-length(z)
    if (length(a)==0) return(list(y=z,l=0,f=0))
    if (is.vector(a)) a<-matrix(a,1,length(a))
    indi<-mkIndi(a,n)
    h<-crossprod(indi,w%*%indi); r<-drop(crossprod(indi,w%*%z))
    b<-solve(h,r); y<-drop(indi%*%b); gy<-2*drop(w%*%(y-z))
    lbd<-mkLagrange(a,gy)
    f<-sum(w*outer(y-z,y-z))
    return(list(y=y,lbd=lbd,f=f,gy=gy))
}