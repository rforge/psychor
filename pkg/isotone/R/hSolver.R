# Huber Loss

hSolver<-function(x,a,extra) {
    w<-extra$w; z<-extra$z; eps<-extra$eps
    fobj<-function(x) sum(w*ifelse(abs(x-z)<2*eps,((x-z)^2)/(4*eps),abs(x-z)-eps))
    gobj<-function(x) w*ifelse(abs(x-z)<2*eps,((x-z))/(2*eps),sign(x-z))
    return(fSolver(x,a,list(fobj=fobj,gobj=gobj)))
}