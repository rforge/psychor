`weighted.pth.fractile` <-
function(y, w)
{
  a <- 1        #alternatively as arguments
  b <- 1
  ox <- order(y)
  y <- y[ox]
  w <- w[ox]
  k <- 1
  low <- cumsum(c(0,w))
  up <- sum(w)-low
  df <- a*low-b*up
  repeat{
	 if (df[k] < 0) k<-k+1
  		else if (df[k] == 0) return((w[k]*y[k]+w[k-1]*y[k-1])/(w[k]+w[k-1]))
  			else return(y[k-1])
	}
}

