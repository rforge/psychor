# Poisson Likelihood

sSolver<-function(z, a, extra) {
    x <- z
    z<-extra$y
    fobj<-function(x) sum(x-z*log(x))
    gobj<-function(x) 1-z/x
    return(fSolver(x,a,list(fobj=fobj,gobj=gobj)))
}