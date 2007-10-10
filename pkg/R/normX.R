`normX` <-
function(x,w) {
  qq<-La.svd(sqrt(w)*x); list(q=(1/sqrt(w))*(qq$u),r=qq$d)
}

