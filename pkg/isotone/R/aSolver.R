# Asymmetric Least Squares

aSolver<-function(x,a,extra) {
    w <- extra$w 
    z <- extra$z
    aw <- extra$aw
    bw <- extra$bw
    fobj <- function(x) sum(w*(x-z)^2*ifelse(x<z, aw, bw))
    gobj <- function(x) 2*w*(x-z)*ifelse(x<z, aw, bw)
    return(fSolver(x,a,list(fobj=fobj,gobj=gobj)))
}