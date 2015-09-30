akerd<-function(xx,hval=NA,aval=.5,op=1,fr=.8,pyhat=FALSE,pts=NA,plotit=FALSE,
xlab="",ylab="",zlab="",theta=50,phi=25,expand=.5,scale=TRUE,ticktype="simple",color='black'){
#
# Compute adaptive kernel density estimate
# 
# (See Silverman, 1986)
#
# op=1 Use expected frequency as initial estimate of the density
# op=2 Univariate case only
#      Use normal kernel to get initial estimate of the density
#  ticktype="detailed" will create ticks as done for a two-dimensional plot
#
#  Note, when pyhat=T, returns estimate of density at points if pts AFTER
#  putting the points in ascending order.
#
xx=elimna(xx)
fval<-"Done"
if(is.matrix(xx)){
if(ncol(xx)>1)fval<-akerdmul(xx,pts=pts,hval=hval,aval=aval,fr=fr,pr=pyhat,
plotit=plotit,theta=theta,phi=phi,expand=expand,scale=scale,ticktype=ticktype)
plotit<-F
}
if(is.matrix(xx) && ncol(xx)==1)xx<-xx[,1]
if(!is.matrix(xx)){
x<-sort(xx)
if(op==1){
m<-mad(x)
if(m==0){
temp<-idealf(x)
m<-(temp$qu-temp$ql)/(qnorm(.75)-qnorm(.25))
}
if(m==0)m<-sqrt(winvar(x)/.4129)
if(m==0)stop("All measures of dispersion are equal to 0")
fhat <- rdplot(x,pyhat=TRUE,plotit=FALSE,fr=fr)
if(m>0)fhat<-fhat/(2*fr*m)
}
if(op==2){
init<-density(xx)
fhat <- init$y
x<-init$x
}
n<-length(x)
if(is.na(hval)){
sig<-sqrt(var(x))
temp<-idealf(x)
iqr<-(temp$qu-temp$ql)/1.34
A<-min(c(sig,iqr))
if(A==0)A<-sqrt(winvar(x))/.64
hval<-1.06*A/length(x)^(.2)
# See Silverman, 1986, pp. 47-48
}
gm<-exp(mean(log(fhat[fhat>0])))
alam<-(fhat/gm)^(0-aval)
dhat<-NA
if(is.na(pts[1]))pts<-x
pts<-sort(pts)
for(j in 1:length(pts)){
temp<-(pts[j]-x)/(hval*alam)
epan<-ifelse(abs(temp)<sqrt(5),.75*(1-.2*temp^2)/sqrt(5),0)
dhat[j]<-mean(epan/(alam*hval))
}
if(plotit){
plot(pts,dhat,type="n",ylab=ylab,xlab=xlab)
lines(pts,dhat,col=color)
}
if(pyhat)fval<-dhat
}
fval
}
