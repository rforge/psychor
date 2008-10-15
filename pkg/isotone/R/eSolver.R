# Approximate l_1

eSolver<-function(x,a,extra) {
    w<-extra$w; z<-extra$z; eps<-extra$eps
    fobj<-function(x) sum(w*sqrt((x-z)^2+eps))
    gobj<-function(x) w*(x-z)/sqrt((x-z)^2+eps)
    return(fSolver(x,a,list(fobj=fobj,gobj=gobj)))
}