# Power Norms

oSolver<-function(x,a,extra) {
    w<-extra$w; z<-extra$z; pow<-extra$p
    fobj<-function(x) sum(w*(abs(x-z)^pow))
    gobj<-function(x) pow*w*sign(x-z)*abs(x-z)^(pow-1)
    return(fSolver(x,a,list(fobj=fobj,gobj=gobj)))
}