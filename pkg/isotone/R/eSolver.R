# Approximate l_1

eSolver<-function(z, a, extra) {
    x <- z
    w <- extra$weights
    z <- extra$y
    eps <- extra$eps
    fobj<-function(x) sum(w*sqrt((x-z)^2+eps))
    gobj<-function(x) w*(x-z)/sqrt((x-z)^2+eps)
    return(fSolver(x,a,list(fobj=fobj,gobj=gobj)))
}