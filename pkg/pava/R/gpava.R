`gpava` <-
function(x, y, w = NULL, solver = weighted.mean, merger = c, ties = "none")
{
# y ... response; either a single vector or a list of vectors (blocks)
# x ... predictor (1 predictor only so far, maybe extension to >1, i.e. generalized pava)
# w ... weights; either a single vector or a list if vectors (weights)
# solver ... either weighted.mean, weighted.median, weighted.pth.fractile, or a user-specified function
# ties ... string for tie treatment: either "primary", "secondary", "tertiary".
    
    
    n <- length(y)
    if(is.null(w)) {
        w <- if(is.list(y))
            lapply(sapply(y, length), function(u) rep.int(1, u))
        else
            rep.int(1, n)
    } else if(is.list(y)) 
        w <- as.list(w)
    y1 <- y
    
 #------------- ties -----------------
 #!!!FIXME: so far this works for y-vectors only!!! to be extended to lists
 if ((is.list(y)) && (ties != "none")) ties <- "none" 
  
    if (ties == "primary") {
      o <- order(x,y)
      r <- order(o)
      y <- y[o]
      w <- w[o]
    }
    if (ties == "secondary") {
      wag <- tapply(w,x,sum) 
      yag <- tapply(y,x,mean)
      xag <- tapply(x,x,mean) 
      o <- order(xag)
      r <- order(o)
      y <- yag[o]
      w <- wag[o]
    }
    if (ties == "tertiary") {
      wag <- tapply(w,x,sum)
      yag <- tapply(y,x,mean)
      xag <- tapply(x,x,mean)
      o <- order(xag)
      r <- order(o)
      y <- yag[o]
      w <- wag[o]  
    }
  #----------- end ties --------------    
    n <- length(y)
    if(is.null(w)) {
        w <- if(is.list(y))
            lapply(sapply(y, length), function(u) rep.int(1, u))
        else
            rep.int(1, n)
    } else if(is.list(y)) 
        w <- as.list(w)
    inds <- as.list(seq_len(n))    
    vals <- mapply(solver, y, w)
    

    ## Combine blocks i and i + 1.    
    combine <- if(is.list(y)) {
        ## In the repeated data case, we explicitly merge the data (and
        ## weight) lists.
        function(i) {
            ## Merge the data and indices, solve, and put things back
            ## into position i, dropping position i + 1.
            j <- i + 1L
            y[[i]] <<- merger(y[[i]], y[[j]])
            w[[i]] <<- c(w[[i]], w[[j]])
            vals[i] <<- solver(y[[i]], w[[i]])
            inds[[i]] <<- c(inds[[i]], inds[[j]])
            keep <- seq_len(n)[-j]
            y <<- y[keep]
            w <<- w[keep]
            vals <<- vals[keep]
            inds <<- inds[keep]
            n <<- n - 1L
        }
    } else {
        function(i) {
            ## In the "simple" case, merge only indices and values.
            j <- i + 1L
            inds[[i]] <<- c(inds[[i]], inds[[j]])
            vals[i] <<- solver(y[inds[[i]]], w[inds[[i]]])
            keep <- seq_len(n)[-j]
            vals <<- vals[keep]
            inds <<- inds[keep]
            n <<- n - 1L
        }
    }
        
    i <- 1L
    repeat {
        if(i < n) {
            if((vals[i] > vals[i + 1])) {
                combine(i)
                while((i > 1L) && (vals[i - 1L] > vals[i])) {
                    combine(i - 1L)
                    i <- i - 1L
                }
            }
            else
                i <- i + 1L
            }
        else break
    }
 yfit.notie <- rep.int(vals, sapply(inds, length))
 
 if (ties == "none") yfit <- yfit.notie
 if (ties == "primary") yfit <- yfit.notie[r]
 if (ties == "secondary") yfit <- as.vector(ifelse(outer(x,xag,"=="),1,0)%*%yfit.notie[r])
 if (ties == "tertiary") yfit <- as.vector(y1 + ifelse(outer(x,xag,"=="),1,0)%*%(yfit.notie[r]-yag[o]))
 
 result <- list(yfit = yfit, y = y1, x = x)  
 class(result) <- "pava"
 result

}

