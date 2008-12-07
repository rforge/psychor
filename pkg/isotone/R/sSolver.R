# Poisson Likelihood

sSolver <- function(z, a, extra) {
    x <- z
    z <- extra$y
    fobj <- function(x) sum(x-z*log(x))       #target value
    gobj <- function(x) 1-z/x                 #gradient
    return(fSolver(x,a,list(fobj=fobj,gobj=gobj)))
}