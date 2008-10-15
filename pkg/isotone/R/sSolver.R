# Poisson Likelihood

sSolver<-function(x,a,extra) {
    z<-extra$z
    fobj<-function(x) sum(x-z*log(x))
    gobj<-function(x) 1-z/x
    return(fSolver(x,a,list(fobj=fobj,gobj=gobj)))
}