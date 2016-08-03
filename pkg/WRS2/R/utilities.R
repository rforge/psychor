matl<-function(x){
  #
  # take data in list mode and store it in a matrix
  #
  J=length(x)
  nval=NA
  for(j in 1:J)nval[j]=length(x[[j]])
  temp<-matrix(NA,ncol=J,nrow=max(nval))
  for(j in 1:J)temp[1:nval[j],j]<-x[[j]]
  temp
}

list2mat=matl

list2vec<-function(x){
  if(!is.list(x))stop("x should have list mode")
  res=as.vector(matl(x))
  res
}

discANOVA.sub<-function(x){
  #
  #
  x=lapply(x,elimna)
  vals=lapply(x,unique)
  vals=sort(elimna(unique(list2vec(vals))))
  n=lapply(x,length)
  n=list2vec(n)
  K=length(vals)
  J=length(x)
  C1=matrix(0,nrow=K,ncol=J)
  for(j in 1:J){
    for(i in 1:K){
      C1[i,j]=C1[i,j]+sum(x[[j]]==vals[i])
    }
    C1[,j]=C1[,j]/n[j]
  }
  test=0
  for(i in 1:K)test=test+var(C1[i,])
  list(test=test,C1=C1)
}

modgen<-function(p,adz=FALSE){
  #
  # Used by regpre to generate all models
  # p=number of predictors
  # adz=T, will add the model where only a measure
  # of location is used.
  #
  #
  model<-list()
  if(p>5)stop("Current version is limited to 5 predictors")
  if(p==1)model[[1]]<-1
  if(p==2){
    model[[1]]<-1
    model[[2]]<-2
    model[[3]]<-c(1,2)
  }
  if(p==3){
    for(i in 1:3)model[[i]]<-i
    model[[4]]<-c(1,2)
    model[[5]]<-c(1,3)
    model[[6]]<-c(2,3)
    model[[7]]<-c(1,2,3)
  }
  if(p==4){
    for(i in 1:4)model[[i]]<-i
    model[[5]]<-c(1,2)
    model[[6]]<-c(1,3)
    model[[7]]<-c(1,4)
    model[[8]]<-c(2,3)
    model[[9]]<-c(2,4)
    model[[10]]<-c(3,4)
    model[[11]]<-c(1,2,3)
    model[[12]]<-c(1,2,4)
    model[[13]]<-c(1,3,4)
    model[[14]]<-c(2,3,4)
    model[[15]]<-c(1,2,3,4)
  }
  if(p==5){
    for(i in 1:5)model[[i]]<-i
    model[[6]]<-c(1,2)
    model[[7]]<-c(1,3)
    model[[8]]<-c(1,4)
    model[[9]]<-c(1,5)
    model[[10]]<-c(2,3)
    model[[11]]<-c(2,4)
    model[[12]]<-c(2,5)
    model[[13]]<-c(3,4)
    model[[14]]<-c(3,5)
    model[[15]]<-c(4,5)
    model[[16]]<-c(1,2,3)
    model[[17]]<-c(1,2,4)
    model[[18]]<-c(1,2,5)
    model[[19]]<-c(1,3,4)
    model[[20]]<-c(1,3,5)
    model[[21]]<-c(1,4,5)
    model[[22]]<-c(2,3,4)
    model[[23]]<-c(2,3,5)
    model[[24]]<-c(2,4,5)
    model[[25]]<-c(3,4,5)
    model[[26]]<-c(1,2,3,4)
    model[[27]]<-c(1,2,3,5)
    model[[28]]<-c(1,2,4,5)
    model[[29]]<-c(1,3,4,5)
    model[[30]]<-c(2,3,4,5)
    model[[31]]<-c(1,2,3,4,5)
  }
  if(adz){
    ic<-length(model)+1
    model[[ic]]<-0
  }
  model
}

pb2gen1<-function(x,y,alpha=.05,nboot=2000,est=onestep,SEED=TRUE,pr=FALSE,...){
  #
  #   Compute a bootstrap confidence interval for the
  #   the difference between any two parameters corresponding to
  #   independent groups.
  #   By default, M-estimators are compared.
  #   Setting est=mean, for example, will result in a percentile
  #   bootstrap confidence interval for the difference between means.
  #   Setting est=onestep will compare M-estimators of location.
  #   The default number of bootstrap samples is nboot=2000
  #
  x<-x[!is.na(x)] # Remove any missing values in x
  y<-y[!is.na(y)] # Remove any missing values in y
  if(SEED)set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  datax<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
  datay<-matrix(sample(y,size=length(y)*nboot,replace=TRUE),nrow=nboot)
  bvecx<-apply(datax,1,est,...)
  bvecy<-apply(datay,1,est,...)
  bvec<-sort(bvecx-bvecy)
  low<-round((alpha/2)*nboot)+1
  up<-nboot-low
  temp<-sum(bvec<0)/nboot+sum(bvec==0)/(2*nboot)
  sig.level<-2*(min(temp,1-temp))
  se<-var(bvec)
  list(est.1=est(x,...),est.2=est(y,...),est.dif=est(x,...)-est(y,...),ci=c(bvec[low],bvec[up]),p.value=sig.level,sq.se=se,n1=length(x),n2=length(y))
}

bootdpci<-function(x,y,est=onestep,nboot=NA,alpha=.05,plotit=TRUE,dif=TRUE,BA=FALSE,SR=TRUE,...){
  #
  #   Use percentile bootstrap method,
  #   compute a .95 confidence interval for the difference between
  #   a measure of location or scale
  #   when comparing two dependent groups.
  #   By default, a one-step M-estimator (with Huber's psi) is used.
  #   If, for example, it is desired to use a fully iterated
  #   M-estimator, use fun=mest when calling this function.
  #
  okay=FALSE
  if(identical(est,onestep))okay=TRUE
  if(identical(est,mom))okay=TRUE
  if(!okay)SR=FALSE
  output<-rmmcppb(x,y,est=est,nboot=nboot,alpha=alpha,SR=SR,
                  plotit=plotit,dif=dif,BA=BA,...)$output
  list(output=output)
}



rmmcppb<-function(x,y=NULL,alpha=.05,
                  con=0,est=onestep,plotit=TRUE,dif=TRUE,grp=NA,nboot=NA,BA=FALSE,hoch=FALSE,xlab="Group 1",ylab="Group 2",pr=TRUE,SEED=TRUE,SR=FALSE,...){
  #
  #   Use a percentile bootstrap method to  compare dependent groups.
  #   By default,
  #   compute a .95 confidence interval for all linear contrasts
  #   specified by con, a J-by-C matrix, where  C is the number of
  #   contrasts to be tested, and the columns of con are the
  #   contrast coefficients.
  #   If con is not specified, all pairwise comparisons are done.
  #
  #   If est=onestep or mom, method SR (see my book on robust methods)
  #   is used to control the probability of at least one Type I error.
  #
  #   Otherwise, Hochberg is used.
  #
  #   dif=T indicates that difference scores are to be used
  #   dif=F indicates that measure of location associated with
  #   marginal distributions are used instead.
  #
  #   nboot is the bootstrap sample size. If not specified, a value will
  #   be chosen depending on the number of contrasts there are.
  #
  #   x can be an n by J matrix or it can have list mode
  #   for two groups, data for second group can be put in y
  #   otherwise, assume x is a matrix (n by J) or has list mode.
  #
  #   A sequentially rejective method is used to control alpha using method SR.
  #
  #   Argument BA: When using dif=F, BA=T uses a correction term
  #   when computing a p-value.
  #
  if(hoch)SR=FALSE #Assume Hochberg if hoch=TRUE even if SR=TRUE
  if(SR){
    okay=FALSE
    if(identical(est,onestep))okay=TRUE
    if(identical(est,mom))okay=TRUE
    SR=okay # 'Only use method SR (argument SR=T) when est=onestep or mom
  }
  if(dif){
    if(pr)print("dif=T, so analysis is done on difference scores")
    temp<-rmmcppbd(x,y=y,alpha=.05,con=con,est,plotit=plotit,grp=grp,nboot=nboot,
                   hoch=TRUE,...)
    output<-temp$output
    con<-temp$con
  }
  if(!dif){
    if(pr){
      print("dif=F, so analysis is done on marginal distributions")
      if(!BA){
        if(identical(est,onestep))print("With M-estimator or MOM, suggest using BA=T and hoch=T")
        if(identical(est,mom))print("With M-estimator or MOM, suggest using BA=T and hoch=T")
      }}
    if(!is.null(y[1]))x<-cbind(x,y)
    if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
    if(is.list(x)){
      if(is.matrix(con)){
        if(length(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
      }}
    if(is.list(x)){
      # put the data in an n by J matrix
      mat<-matl(x)
    }
    if(is.matrix(x) && is.matrix(con)){
      if(ncol(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
      mat<-x
    }
    if(is.matrix(x))mat<-x
    if(!is.na(sum(grp)))mat<-mat[,grp]
    mat<-elimna(mat) # Remove rows with missing values.
    x<-mat
    J<-ncol(mat)
    xcen<-x
    for(j in 1:J)xcen[,j]<-x[,j]-est(x[,j],...)
    Jm<-J-1
    if(sum(con^2)==0){
      d<-(J^2-J)/2
      con<-matrix(0,J,d)
      id<-0
      for (j in 1:Jm){
        jp<-j+1
        for (k in jp:J){
          id<-id+1
          con[j,id]<-1
          con[k,id]<-0-1
        }}}
    d<-ncol(con)
    if(is.na(nboot)){
      if(d<=4)nboot<-1000
      if(d>4)nboot<-5000
    }
    n<-nrow(mat)
    crit.vec<-alpha/c(1:d)
    connum<-ncol(con)
    if(SEED)set.seed(2) # set seed of random number generator so that
    #             results can be duplicated.
    xbars<-apply(mat,2,est,...)
    psidat<-NA
    for (ic in 1:connum)psidat[ic]<-sum(con[,ic]*xbars)
    psihat<-matrix(0,connum,nboot)
    psihatcen<-matrix(0,connum,nboot)
    bvec<-matrix(NA,ncol=J,nrow=nboot)
    bveccen<-matrix(NA,ncol=J,nrow=nboot)
    if(pr)print("Taking bootstrap samples. Please wait.")
    data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
    for(ib in 1:nboot){
      bvec[ib,]<-apply(x[data[ib,],],2,est,...)
      bveccen[ib,]<-apply(xcen[data[ib,],],2,est,...)
    }
    #
    # Now have an nboot by J matrix of bootstrap values.
    #
    test<-1
    bias<-NA
    for (ic in 1:connum){
      psihat[ic,]<-apply(bvec,1,bptdpsi,con[,ic])
      psihatcen[ic,]<-apply(bveccen,1,bptdpsi,con[,ic])
      bias[ic]<-sum((psihatcen[ic,]>0))/nboot-.5
      ptemp<-(sum(psihat[ic,]>0)+.5*sum(psihat[ic,]==0))/nboot
      if(BA)test[ic]<-ptemp-.1*bias[ic]
      if(!BA)test[ic]<-ptemp
      test[ic]<-min(test[ic],1-test[ic])
      test[ic]<-max(test[ic],0)  # bias corrected might be less than zero
    }
    test<-2*test
    ncon<-ncol(con)
    dvec<-alpha/c(1:ncon) # Assume Hochberg unless specified otherwise
    if(SR){
      if(alpha==.05){
        dvec<-c(.025,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
        dvecba<-c(.05,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
        if(ncon > 10){
          avec<-.05/c(11:ncon)
          dvec<-c(dvec,avec)
        }}
      if(alpha==.01){
        dvec<-c(.005,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
        dvecba<-c(.01,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
        if(ncon > 10){
          avec<-.01/c(11:ncon)
          dvec<-c(dvec,avec)
        }}
      if(alpha != .05 && alpha != .01){
        dvec<-alpha/c(1:ncon)
        dvecba<-dvec
        dvec[2]<-alpha
      }}
    if(hoch)dvec<-alpha/c(1:ncon)
    dvecba<-dvec
    if(plotit && ncol(bvec)==2){
      z<-c(0,0)
      one<-c(1,1)
      plot(rbind(bvec,z,one),xlab=xlab,ylab=ylab,type="n")
      points(bvec)
      totv<-apply(x,2,est,...)
      cmat<-var(bvec)
      dis<-mahalanobis(bvec,totv,cmat)
      temp.dis<-order(dis)
      ic<-round((1-alpha)*nboot)
      xx<-bvec[temp.dis[1:ic],]
      xord<-order(xx[,1])
      xx<-xx[xord,]
      temp<-chull(xx)
      lines(xx[temp,])
      lines(xx[c(temp[1],temp[length(temp)]),])
      abline(0,1)
    }
    temp2<-order(0-test)
    ncon<-ncol(con)
    zvec<-dvec[1:ncon]
    if(BA)zvec<-dvecba[1:ncon]
    sigvec<-(test[temp2]>=zvec)
    output<-matrix(0,connum,6)
    dimnames(output)<-list(NULL,c("con.num","psihat","p.value","p.sig","ci.lower","ci.upper"))
    tmeans<-apply(mat,2,est,...)
    psi<-1
    output[temp2,4]<-zvec
    for (ic in 1:ncol(con)){
      output[ic,2]<-sum(con[,ic]*tmeans)
      output[ic,1]<-ic
      output[ic,3]<-test[ic]
      temp<-sort(psihat[ic,])
      #icl<-round(output[ic,4]*nboot/2)+1
      icl<-round(alpha*nboot/2)+1
      icu<-nboot-(icl-1)
      output[ic,5]<-temp[icl]
      output[ic,6]<-temp[icu]
    }
  }
  num.sig<-sum(output[,3]<=output[,4])
  list(output=output,con=con,num.sig=num.sig)
}

bptdpsi<-function(x,con){
  # Used by bptd to compute bootstrap psihat values
  #
  bptdpsi<-sum(con*x)
  bptdpsi
}

rmmismcp<-function(x,y=NA,alpha=.05,con=0,est=tmean,plotit=TRUE,grp=NA,nboot=500,
                   SEED=TRUE,xlab="Group 1",ylab="Group 2",pr=FALSE,...){
  #
  #   Use a percentile bootstrap method to  compare  marginal measures of location for dependent groups.
  #   Missing values are allowed; vectors of observations that contain
  #   missing values are not simply removed as done by rmmcppb.
  #   Only marginal measures of location are compared,
  #   The function computes a .95 confidence interval for all linear contrasts
  #   specified by con, a J by C matrix, where  C is the number of
  #   contrasts to be tested, and the columns of con are the
  #   contrast coefficients.
  #   If con is not specified, all pairwise comparisons are done.
  #
  #   By default, a 20% trimmed is used and a sequentially rejective method
  #   is used to control the probability of at least one Type I error.
  #
  #   nboot is the bootstrap sample size.
  #
  #   x can be an n by J matrix or it can have list mode
  #   for two groups, data for second group can be put in y
  #   otherwise, assume x is a matrix (n by J) or has list mode.
  #
  #
  if(!is.na(y[1]))x<-cbind(x,y)
  if(is.list(x))x=matl(x)
  if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
  if(is.list(x)){
    if(is.matrix(con)){
      if(length(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
    }}
  if(is.list(x)){
    # put the data in an n by J matrix
    mat<-matl(x)
  }
  if(is.matrix(x) && is.matrix(con)){
    if(ncol(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
    mat<-x
  }
  J<-ncol(x)
  Jm<-J-1
  flag.con=F
  if(sum(con^2)==0){
    flag.con=T
    d<-(J^2-J)/2
    con<-matrix(0,J,d)
    id<-0
    for (j in 1:Jm){
      jp<-j+1
      for (k in jp:J){
        id<-id+1
        con[j,id]<-1
        con[k,id]<-0-1
      }}}
  d<-ncol(con)
  n<-nrow(x)
  crit.vec<-alpha/c(1:d)
  connum<-ncol(con)
  if(SEED)set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  xbars<-apply(x,2,est,na.rm=TRUE)
  psidat<-NA
  bveccen<-matrix(NA,ncol=J,nrow=nboot)
  for (ic in 1:connum)psidat[ic]<-sum(con[,ic]*xbars)
  psihat<-matrix(0,connum,nboot)
  psihatcen<-matrix(0,connum,nboot)
  bvec<-matrix(NA,ncol=J,nrow=nboot)
  if(pr)print("Taking bootstrap samples. Please wait.")
  data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
  for(ib in 1:nboot){
    bvec[ib,]<-apply(x[data[ib,],],2,est,na.rm=TRUE,...)
  }
  #
  # Now have an nboot by J matrix of bootstrap measures of location.
  #
  test<-1
  for (ic in 1:connum){
    for(ib in 1:nboot){
      psihat[ic,ib]=sum(con[,ic]*bvec[ib,])
    }
    matcon=c(0,psihat[ic,])
    dis=mean((psihat[ic,]<0))+.5*mean((psihat[ic,]==0))
    test[ic]<-2*min(c(dis,1-dis)) # the p-value
  }
  ncon<-ncol(con)
  dvec<-alpha/c(1:ncon)
  if(plotit && ncol(bvec)==2){
    z<-c(0,0)
    one<-c(1,1)
    plot(rbind(bvec,z,one),xlab=xlab,ylab=ylab,type="n")
    points(bvec)
    totv<-apply(x,2,est,na.rm=TRUE,...)
    cmat<-var(bvec)
    dis<-mahalanobis(bvec,totv,cmat)
    temp.dis<-order(dis)
    ic<-round((1-alpha)*nboot)
    xx<-bvec[temp.dis[1:ic],]
    xord<-order(xx[,1])
    xx<-xx[xord,]
    temp<-chull(xx)
    lines(xx[temp,])
    lines(xx[c(temp[1],temp[length(temp)]),])
    abline(0,1)
  }
  temp2<-order(0-test)
  ncon<-ncol(con)
  zvec<-dvec[1:ncon]
  sigvec<-(test[temp2]>=zvec)
  output<-matrix(0,connum,6)
  dimnames(output)<-list(NULL,c("con.num","psihat","p.value",
                                "crit.sig","ci.lower","ci.upper"))
  tmeans<-apply(x,2,est,na.rm=TRUE,...)
  psi<-1
  output[temp2,4]<-zvec
  for (ic in 1:ncol(con)){
    output[ic,2]<-sum(con[,ic]*tmeans)
    output[ic,1]<-ic
    output[ic,3]<-test[ic]
    temp<-sort(psihat[ic,])
    icl<-round(output[ic,4]*nboot/2)+1
    icu<-nboot-(icl-1)
    output[ic,5]<-temp[icl]
    output[ic,6]<-temp[icu]
  }
  if(!flag.con){
  }
  if(flag.con){
    CC=(J^2-J)/2
    test<-matrix(NA,CC,7)
    dimnames(test)<-list(NULL,c("Group","Group","psi.hat","p.value","p.crit",
                                "ci.low","ci.upper"))
    jcom<-0
    for (j in 1:J){
      for (k in 1:J){
        if (j < k){
          jcom<-jcom+1
          test[jcom,1]=j
          test[jcom,2]=k
          test[jcom,3:5]=output[jcom,2:4]
          test[jcom,6:7]=output[jcom,5:6]
          con=NULL
        }}}}
  if(!flag.con)test=output
  #num.sig<-sum(output[,4]<=output[,5])
  if(flag.con)num.sig<-sum(test[,4]<=test[,5])
  if(!flag.con)num.sig<-sum(test[,3]<=test[,4])
  list(output=test,con=con,num.sig=num.sig)
}


rmmcppbd<-function(x,y=NULL,alpha=.05,con=0,est=onestep,plotit=TRUE,grp=NA,nboot=NA,
                   hoch=TRUE,SEED=TRUE,...){
  #
  #   Use a percentile bootstrap method to  compare dependent groups
  #   based on difference scores.
  #   By default,
  #   compute a .95 confidence interval for all linear contrasts
  #   specified by con, a J by C matrix, where  C is the number of
  #   contrasts to be tested, and the columns of con are the
  #   contrast coefficients.
  #   If con is not specified, all pairwise comparisons are done.
  #
  #   By default, one-step M-estimator is used
  #    and a sequentially rejective method
  #   is used to control the probability of at least one Type I error.
  #
  #   nboot is the bootstrap sample size. If not specified, a value will
  #   be chosen depending on the number of contrasts there are.
  #
  #   x can be an n by J matrix or it can have list mode
  #   for two groups, data for second group can be put in y
  #   otherwise, assume x is a matrix (n by J) or has list mode.
  #
  #   A sequentially rejective method is used to control alpha.
  #   If n>=80, hochberg's method is used.
  #
  if(!is.null(y[1]))x<-cbind(x,y)
  if(!is.list(x) && !is.matrix(x))stop("Data must be stored in a matrix or in list mode.")
  if(is.list(x)){
    if(is.matrix(con)){
      if(length(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
    }}
  if(is.list(x)){
    # put the data in an n by J matrix
    mat<-matl(x)
  }
  if(is.matrix(x) && is.matrix(con)){
    if(ncol(x)!=nrow(con))stop("The number of rows in con is not equal to the number of groups.")
    mat<-x
  }
  if(is.matrix(x))mat<-x
  if(!is.na(sum(grp)))mat<-mat[,grp]
  x<-mat
  mat<-elimna(mat) # Remove rows with missing values.
  x<-mat
  J<-ncol(mat)
  n=nrow(mat)
  if(n>=80)hoch=T
  Jm<-J-1
  if(sum(con^2)==0){
    d<-(J^2-J)/2
    con<-matrix(0,J,d)
    id<-0
    for (j in 1:Jm){
      jp<-j+1
      for (k in jp:J){
        id<-id+1
        con[j,id]<-1
        con[k,id]<-0-1
      }}}
  d<-ncol(con)
  if(is.na(nboot)){
    nboot<-5000
    if(d<=10)nboot<-3000
    if(d<=6)nboot<-2000
    if(d<=4)nboot<-1000
  }
  n<-nrow(mat)
  crit.vec<-alpha/c(1:d)
  connum<-ncol(con)
  # Create set of differences based on contrast coefficients
  xx<-x%*%con
  xx<-as.matrix(xx)
  if(SEED)set.seed(2) # set seed of random number generator so that
  #             results can be duplicated.
  psihat<-matrix(0,connum,nboot)
  bvec<-matrix(NA,ncol=connum,nrow=nboot)
  data<-matrix(sample(n,size=n*nboot,replace=TRUE),nrow=nboot)
  # data is an nboot by n matrix
  if(ncol(xx)==1){
    for(ib in 1:nboot)psihat[1,ib]<-est(xx[data[ib,]],...)
  }
  if(ncol(xx)>1){
    for(ib in 1:nboot)psihat[,ib]<-apply(xx[data[ib,],],2,est,...)
  }
  #
  # Now have an nboot by connum matrix of bootstrap values.
  #
  test<-1
  for (ic in 1:connum){
    test[ic]<-(sum(psihat[ic,]>0)+.5*sum(psihat[ic,]==0))/nboot
    test[ic]<-min(test[ic],1-test[ic])
  }
  test<-2*test
  ncon<-ncol(con)
  if(alpha==.05){
    dvec<-c(.025,.025,.0169,.0127,.0102,.00851,.0073,.00639,.00568,.00511)
    if(ncon > 10){
      avec<-.05/c(11:ncon)
      dvec<-c(dvec,avec)
    }}
  if(alpha==.01){
    dvec<-c(.005,.005,.00334,.00251,.00201,.00167,.00143,.00126,.00112,.00101)
    if(ncon > 10){
      avec<-.01/c(11:ncon)
      dvec<-c(dvec,avec)
    }}
  if(alpha != .05 && alpha != .01){
    dvec<-alpha/c(1:ncon)
    dvec[2]<-alpha/2
  }
  if(hoch)dvec<-alpha/(2*c(1:ncon))
  dvec<-2*dvec
  if(plotit && connum==1){
    plot(c(psihat[1,],0),xlab="",ylab="Est. Difference")
    points(psihat[1,])
    abline(0,0)
  }
  temp2<-order(0-test)
  ncon<-ncol(con)
  zvec<-dvec[1:ncon]
  sigvec<-(test[temp2]>=zvec)
  output<-matrix(0,connum,6)
  dimnames(output)<-list(NULL,c("con.num","psihat","p.value","p.crit","ci.lower","ci.upper"))
  tmeans<-apply(xx,2,est,...)
  psi<-1
  icl<-round(dvec[ncon]*nboot/2)+1
  icu<-nboot-icl-1
  for (ic in 1:ncol(con)){
    output[ic,2]<-tmeans[ic]
    output[ic,1]<-ic
    output[ic,3]<-test[ic]
    output[temp2,4]<-zvec
    temp<-sort(psihat[ic,])
    output[ic,5]<-temp[icl]
    output[ic,6]<-temp[icu]
  }
  num.sig<-sum(output[,3]<=output[,4])
  list(output=output,con=con,num.sig=num.sig)
}




