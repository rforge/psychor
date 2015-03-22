FunDispariSEM3 <-
function(Xi, Xim, R, DXm, SSk, theta0, saturado) {
## Input parameters:
    ## Xi matrix of effect indicators as column vectors
    ## Xim centered effect indicators
    ## R number of matrices
    ## DXm centered distances in column vector
    ## SSk decreasing constraints for identification 
    ## theta0 initial parameter values in SEM
    ## saturado indicator value of satured model (0 is the default not satured)
    
## Output values:
    ## Delta Least squares estimated symmetric disparity values
    ## theta vector of estimated parameter values in SEM

cotainf <- rep(0, R+3)
cotasup <- rep(Inf, R+3)

theta <- NULL         
if ((saturado == 1) && (R == 2)) { ## Satured model (only for R=2)
    S <- cov(cbind(DXm, Xim))
    x1 <- Xi[,1]
    x2 <- Xi[,2]
    theta[1] <- sqrt(S[2,1]*S[3,1]/S[3,2])
    theta[2] <- sqrt(S[3,2]*S[2,1]/S[3,1])
    theta[3] <- sqrt(S[3,2]*S[3,1]/S[2,1])
    theta[4] <- S[2,2]-theta[2]^2
    theta[5] <- S[3,3]-theta[3]^2
    ## FIXME: check saturated
    Delta <- c(theta[2], theta[3]) %*% solve(rbind(c(theta[2]^2+theta[4], theta[2]*theta[3]), c(theta[2]*theta[3], theta[3]^2+theta[5]))) %*% rbind(x1, x2)
    Delta <- t(Delta)
} else  {                                     ## Not satured model
    
    ## FIXME!!! 
    ## optimset('Display','Final','TolFun',0.0000001,'MaxIter',1000);
    ## [theta] = lsqnonlin(@(theta)ULS3(theta,Xim,DXm,SSk),theta0,cotainf,cotasup); ## Row vector.
    ## see nonlin_residmin (@ (p) ..., x0, optimset ("lbound", lb, "ubound", ub))
    theta <- lsqnonlin(ULS3, theta0, Xim = Xim, DXm = DXm, SSk = SSk, options = list(maxeval = 1000, tolg = 1e-7))$x
     
    Lambda <- theta[2:(R+1)]  ## Row vector.
    Sgm <- matrix(0, R, R)
    for (i in (2:(R+1))) Sgm[i-1,i-1] <- theta[i]^2+theta[R+2]
  
    for (i in 3:(R+1)) {
        for (j in 2:(i-1)) {
            Sgm[i-1,j-1] <- theta[i]*theta[j]
            Sgm[j-1,i-1] <- Sgm[i-1,j-1]
        }
    }
    Delta <- Lambda %*% solve(Sgm) %*% t(Xi)
    Delta <- t(Delta)     ## Column vector.
}

result <- list(Delta = Delta, theta = theta)
return(result)
}
