`sumSet` <-
function(g,n,p,y,set,active) {
  z<-array(0.0,c(n,p))
  for (j in set)
  	if (active[j]){
  	    gg<-g[,j]; ii<-which(!is.na(gg));
  		z[ii,]<-z[ii,]+y[[j]][gg[ii],]}
  return(z)
}

