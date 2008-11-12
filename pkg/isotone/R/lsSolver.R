#Least squares with diagonal weights
#returns fitted values y, Lagrange multiplier (lambda) lbd, function value f, gradient gy

lsSolver<-function(x,a,extra) 
{
#x ... function values
#a ... active constraints

    w <- extra$w                                 #weights
    z <- extra$z                                 #response
    n <- length(z)
    if (length(a)==0) return(list(y=z,l=0,f=0))  #no active set, break
    if (is.vector(a)) a <- matrix(a,1,length(a)) #only 1 active set (as matrix)
    
    indi <- mkIndi(a,n)                          #compute indicators
    h <- crossprod(indi, w*indi)
    r <- drop(crossprod(indi,w*z))
    b <- solve(h,r)
    y <- drop(indi%*%b)
    gy <- 2*w*(y-z)                              #gradient
    lbd <- mkLagrange(a, gy)
    f <- sum(w*(y-z)^2)                          #value target function
    return(list(y = y, lbd = lbd, f = f, gy = gy))
}