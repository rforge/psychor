ifmest<-function(x,bend=1.28,op=2){
#
#   Estimate the influence function of an M-estimator, using
#   Huber's Psi, evaluated at x.
#
#   Data are in the vector x, bend is the percentage bend
#
#  op=2, use adaptive kernel estimator
#  otherwise use Rosenblatt's shifted histogram
#
tt<-mest(x,bend)  # Store M-estimate in tt
s<-mad(x)*qnorm(.75)
if(op==2){
val<-akerd(x,pts=tt,plotit=FALSE,pyhat=T)
val1<-akerd(x,pts=tt-s,plotit=FALSE,pyhat=T)
val2<-akerd(x,pts=tt+s,plotit=FALSE,pyhat=T)
}
if(op!=2){
val<-kerden(x,0,tt)
val1<-kerden(x,0,tt-s)
val2<-kerden(x,0,tt+s)
}
ifmad<-sign(abs(x-tt)-s)-(val2-val1)*sign(x-tt)/val
ifmad<-ifmad/(2*.6745*(val2+val1))
y<-(x-tt)/mad(x)
n<-length(x)
b<-sum(y[abs(y)<=bend])/n
a<-hpsi(y,bend)*mad(x)-ifmad*b
ifmest<-a/(length(y[abs(y)<=bend])/n)
ifmest
}