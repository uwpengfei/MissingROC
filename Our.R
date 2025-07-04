meleobj<-function(para, y, R, x, v)             
{  
  xmat=as.matrix(cbind(x,v))
  n=length(y)
  dismat=xmat ##Covariates for disease model
  dismat1=cbind(rep(1,n),dismat)
  
  mismat=xmat ## Covariates for missing model/verification model
  mismat1=cbind(rep(1,n),mismat)
  
  k1=ncol(dismat1)
  k2=ncol(mismat1)
  
  mu =  para[1:k1]
  psi= para[(k1+1):(k1+k2)]
  beta=   para[k1+k2+1]
  
  
  
  p1xv=1-plogis( as.numeric(dismat1%*%mu) )
  
  part1 = sum( R*( y*log( p1xv+1e-30 ) +(1-y)*log(1-p1xv+1e-30) )  )
  
  
  cxv=log(1-p1xv+p1xv*exp(beta)+1e-30 )
  
  txv=as.numeric(mismat1%*%psi)+cxv
  
  part2=sum((1-R)*txv)
  
  
  part3 =-sum( log(1+ exp(txv) ) ) 
  
  out= part1+part2+part3
  -out
}


###Calcualte MLE###
estpara<- function(y, R, x, v)
{  
  
  n=length(y)
  xmat=cbind(x,v)
  xmat=as.matrix(xmat)
  dismat=xmat ##Covariates for disease model
  dismat1=cbind(rep(1,n),dismat)
  
  mismat=xmat ## Covariates for missing model/verification model
  mismat1=cbind(rep(1,n),mismat)
  
  k1=ncol(dismat1)
  k2=ncol(mismat1)
  
  
  mu0=-glm(y[R>0]~dismat[R>0,],family="binomial")$coef
  mu0=as.numeric(mu0)
  phi0=c(0,rep(0,k2)) 
  
  output=c()  
  for(i in 1:16)
  {
    par0=phi0+runif(k2+1,-0.5,0.5)
    para0=c(mu0,par0)
    
    out=nlminb(para0, meleobj, lower= rep(-10, k1+k2+1), upper=rep(10, k1+k2+1),y=y,R=R,x=x,v=v)
    output=rbind(output,c(out$par,out$objective))
  } 
  
  
  para=output[( (1:nrow(output))[output[,k1+k2+2]<=min(output[,k1+k2+2])])[1],]
  mu =  para[1:k1]
  phi= para[(k1+1):(k1+k2+1)]
list(mu=mu,phi=phi)
}


##Calculate the variance estimate for MLE###

varpara=function(y,R,x,v,mu,phi)
{
  n=length(y)
  cova = cbind(x, v)
  cova=as.matrix(cova)
  cova1 = as.matrix(cova)
  cova2 = as.matrix(cova)
  psi=phi[1:( length(phi)-1)]
  beta=phi[length(phi)]	
  mismat1=cbind(rep(1,n),cova2)
  
  # calculate P(Y=1|X,V,R=1)
  w=cbind(rep(1,n),cova1)
  p1xv=1/(1+exp(as.numeric(w%*%mu)))
  
  # calculate g(x,v,r; theta)
  gmat=cbind(rep(1,n),cova1,(R-1))
  mub=c(mu,beta)
  gxvr=1/(1 + exp(as.numeric(gmat%*%mub)))
  
  # calculate g'(x,v,r; theta)
  gxvr1=matrix(rep(-gxvr*(1-gxvr),ncol(gmat)),nrow=n)*gmat
  
  # c(x,v) and t(x,v)
  cxv=log(1-p1xv+p1xv*exp(beta))
  txv=as.numeric(mismat1%*%psi)+cxv
  
  
  # calculate lambda
  lam=mean(gxvr)
  #lam=mean(R)
  
  # calculate score
  smu1=-matrix(rep(R*(y-p1xv),ncol(w)),nrow=n)*w
  smu2=matrix(rep((R-1/(1+exp(txv)))*(p1xv*(1-p1xv))*(exp(beta)-1)/exp(cxv),ncol(w)),nrow=n)*w
  smu=smu1+smu2
  
  sphi=-matrix(rep(R-1/(1+exp(txv)),ncol(w)),nrow=n)*w
  
  sbeta=-(R-1/(1+exp(txv)))*p1xv*exp(beta)/(1-p1xv+p1xv*exp(beta))
  
  # calculate J^-1(theta)
  Jtheta=cov(cbind(smu,sbeta,sphi))
  
  
  # modify the matrix to prevent too little eigenvalues
  eig=eigen(Jtheta)
  ev=eig$values
  evc=eig$vectors
  ev.new=diag(max(ev)/10^7,ncol(Jtheta))+diag(ev,ncol(Jtheta))
  Jtheta=evc%*%ev.new%*%t(evc)
  J_inv=solve(Jtheta)
  J_inv/n
}



##Calculate the proposed AUC estimate and the se##

estauc=function(y,R,x,v,mu,phi)
  {
  n=length(y)
  cova = cbind(x, v)
  cova=as.matrix(cova)
  cova1 = as.matrix(cova)
  cova2 = as.matrix(cova)
  psi=phi[1:( length(phi)-1)]
  beta=phi[length(phi)]	
  mismat1=cbind(rep(1,n),cova2)
  
  # calculate P(Y=1|X,V,R=1)
  w=cbind(rep(1,n),cova1)
  p1xv=1/(1+exp(as.numeric(w%*%mu)))
  
  # calculate g(x,v,r; theta)
  gmat=cbind(rep(1,n),cova1,(R-1))
  mub=c(mu,beta)
  gxvr=1/(1 + exp(as.numeric(gmat%*%mub)))
  
  w1=gxvr/sum(gxvr)
  w0=(1-gxvr)/(sum(1-gxvr))
  
  ##calculate the point estimate or auc##
  A=newauc(x,w0,w1)
  
  
  # calculate g'(x,v,r; theta)
  gxvr1=matrix(rep(-gxvr*(1-gxvr),ncol(gmat)),nrow=n)*gmat
  
  # c(x,v) and t(x,v)
  cxv=log(1-p1xv+p1xv*exp(beta))
  txv=as.numeric(mismat1%*%psi)+cxv
  
  # calculate F0 and F1
  ux=sort(unique(x))
  hatF0=stepfun(ux,c(0,ewcdf(x,w0)))
  hatF1=stepfun(ux,c(0,ewcdf(x,w1)))
  F0=hatF0(x)
  F1=hatF1(x)
  
  # calculate E1,E2,E3
  #e1=apply(gxvr1,MARGIN = 2,FUN = mean)
  gF0=gxvr1*matrix(rep(F0-A,ncol(gmat)),nrow=n)
  gF1=-gxvr1*matrix(rep(1-F1-A,ncol(gmat)),nrow=n)
  e1=apply(gF0,MARGIN = 2,FUN = mean)
  e2=apply(gF1,MARGIN = 2,FUN = mean)
  
  # calculate lambda
  lam=mean(gxvr)
  #lam=mean(R)
  
  # calculate score
  smu1=-matrix(rep(R*(y-p1xv),ncol(w)),nrow=n)*w
  smu2=matrix(rep((R-1/(1+exp(txv)))*(p1xv*(1-p1xv))*(exp(beta)-1)/exp(cxv),ncol(w)),nrow=n)*w
  smu=smu1+smu2
  
  sphi=-matrix(rep(R-1/(1+exp(txv)),ncol(w)),nrow=n)*w
  
  sbeta=-(R-1/(1+exp(txv)))*p1xv*exp(beta)/(1-p1xv+p1xv*exp(beta))
  
  # calculate J^-1(theta)
  Jtheta=cov(cbind(smu,sbeta,sphi))
  
  
  # modify the matrix to prevent too little eigenvalues
  eig=eigen(Jtheta)
  ev=eig$values
  evc=eig$vectors
  ev.new=diag(max(ev)/10^7,ncol(Jtheta))+diag(ev,ncol(Jtheta))
  Jtheta=evc%*%ev.new%*%t(evc)
  J_inv=solve(Jtheta)
  
  # calculate Sigma_Z 
  z1=cbind(smu,sbeta,sphi)
  z2=gxvr*(F0-A)
  z3=(1-gxvr)*(1-F1-A)
  z=cbind(z1,z2,z3)
  Sigma_z=cov(z)
  
  # calculate delta
  k1=1+length(mu)
  k2=length(psi)
  Ik1k2=cbind(diag(k1),matrix(rep(0,k1*k2),nrow=k1))
  delta1=(t(e1)/lam+t(e2)/(1-lam))%*%Ik1k2%*%J_inv
  delta2=1/lam
  delta3=1/(1-lam)
  delta=c(delta1,delta2,delta3)
  
  # variance
 varauc= t(delta)%*%Sigma_z%*%delta/n
 varauc=as.numeric(varauc)
 seauc=sqrt(varauc)
 
c(A,seauc)
 }

###Calculate the proposed ROC(s) estimate and the se at the given s##


estroc=function(s,y,R,x,v,mu,phi)
  {
  n=length(y)
  cova = cbind(x, v)
  cova=as.matrix(cova)
  cova1 = as.matrix(cova)
  cova2 = as.matrix(cova)
  psi=phi[1:( length(phi)-1)]
  beta=phi[length(phi)]	
  mismat1=cbind(rep(1,n),cova2)
  
  # calculate P(Y=1|X,V,R=1)
  w=cbind(rep(1,n),cova1)
  p1xv=1/(1+exp(as.numeric(w%*%mu)))
  
  # calculate g(x,v,r; theta)
  gmat=cbind(rep(1,n),cova1,(R-1))
  mub=c(mu,beta)
  gxvr=1/(1 + exp(as.numeric(gmat%*%mub)))
  w1=gxvr/sum(gxvr)
  w0=(1-gxvr)/(sum(1-gxvr))

  # calculate g'(x,v,r; theta)
  gxvr1=matrix(rep(-gxvr*(1-gxvr),ncol(gmat)),nrow=n)*gmat
  
  # c(x,v) and t(x,v)
  cxv=log(1-p1xv+p1xv*exp(beta))
  txv=as.numeric(mismat1%*%psi)+cxv
  
  # calculate ROC(s), F0 and F1
  ux=sort(unique(x))
  hatF0=stepfun(ux,c(0,ewcdf(x,w0)))
  hatF1=stepfun(ux,c(0,ewcdf(x,w1)))
  rocs=disroc(s,hatF0(ux),hatF1(ux))
  
  F0=hatF0(x)
  F1=hatF1(x)
  
  
  
  # find 1-s quantile for F0
  x_ord=sort(x)
  xi_s=x_ord[which(sort(F0)>=(1-s))[1]]
  F0s=hatF0(xi_s)
  F1s=hatF1(xi_s)
  
  
  # calculate bandwidth
  w0_fi=w0
  w1_fi=w1
  q25_0=x_ord[which(sort(F0)>=0.25)[1]]
  q75_0=x_ord[which(sort(F0)>=0.75)[1]]
  r_0=q75_0-q25_0
  
  q25_1=x_ord[which(sort(F1)>=0.25)[1]]
  q75_1=x_ord[which(sort(F1)>=0.75)[1]]
  r_1=q75_1-q25_1
  
  sig0=sum(x^2*w0_fi/sum(w0_fi))-sum(x*w0_fi/sum(w0_fi))^2
  sig1=sum(x^2*w1_fi/sum(w1_fi))-sum(x*w1_fi/sum(w1_fi))^2
  
  b0=1.06/n^(1/5)*min(sqrt(sig0),r_0/1.34)
  b1=1.06/n^(1/5)*min(sqrt(sig1),r_1/1.34)
  
  f0s=sum(dnorm((xi_s-x)/b0)*w0_fi/sum(w0_fi)/b0)
  f1s=sum(dnorm((xi_s-x)/b1)*w1_fi/sum(w1_fi)/b1)
  
  
  
  # calculate E3,E4
  gF0=gxvr1*matrix(rep((x<=xi_s)-F0s,ncol(gmat)),nrow=n)
  gF1=gxvr1*matrix(rep((x<=xi_s)-F1s,ncol(gmat)),nrow=n)
  e3=apply(gF0,MARGIN = 2,FUN = mean)
  e4=apply(gF1,MARGIN = 2,FUN = mean)
  
  # calculate lambda
  lam=mean(gxvr)
  #lam=mean(R)
  
  # calculate score
  smu1=-matrix(rep(R*(y-p1xv),ncol(w)),nrow=n)*w
  smu2=matrix(rep((R-1/(1+exp(txv)))*(p1xv*(1-p1xv))*(exp(beta)-1)/exp(cxv),ncol(w)),nrow=n)*w
  smu=smu1+smu2
  
  sphi=-matrix(rep(R-1/(1+exp(txv)),ncol(w)),nrow=n)*w
  
  sbeta=-(R-1/(1+exp(txv)))*p1xv*exp(beta)/(1-p1xv+p1xv*exp(beta))
  
  # calculate J^-1(theta)
  Jtheta=cov(cbind(smu,sbeta,sphi))
  
  
  # modify the matrix to prevent too little eigenvalues
  eig=eigen(Jtheta)
  ev=eig$values
  evc=eig$vectors
  ev.new=diag(max(ev)/10^7,ncol(Jtheta))+diag(ev,ncol(Jtheta))
  Jtheta=evc%*%ev.new%*%t(evc)
  J_inv=solve(Jtheta)
  
  
  
  # calculate Sigma_Z 
  z1=cbind(smu,sbeta,sphi)
  z2=(1-gxvr)*((x<=xi_s)-F0s)
  z3=gxvr*((x<=xi_s)-F1s)
  z=cbind(z1,z2,z3)
  Sigma_z=cov(z)
  
  # calculate delta
  k1=1+length(mu)
  k2=length(psi)
  Ik1k2=cbind(diag(k1),matrix(rep(0,k1*k2),nrow=k1))
  delta1=(-t(e4)/lam-t(e3)/(1-lam)*f1s/f0s)%*%Ik1k2%*%J_inv
  delta2=1/(1-lam)*f1s/f0s
  delta3=-1/lam
  delta=c(delta1,delta2,delta3)
  
  # variance
  varroc=t(delta)%*%Sigma_z%*%delta/n
serocs=sqrt(as.numeric(varroc))
c(rocs,serocs)
}





##Our proposed method for estimating AUC and ROC curve at ss##

ourest<- function(y, R, x, v,ss)
{  
  
 
  para=estpara(y, R, x, v)
  mu=as.numeric(para$mu)
  phi=as.numeric(para$phi)
  n=length(y)
  beta=phi[length(phi)]	
  
  gmat=cbind(rep(1,n),cbind(x,v),(R-1))
  gmat=as.matrix(gmat)
  mub=c(mu,beta)
  gxvr=1/(1 + exp(as.numeric(gmat%*%mub)))
  w1=gxvr/sum(gxvr)
  w0=(1-gxvr)/(sum(1-gxvr))
  
 c(newauc(x,w0,w1),newroc(ss,x,w0,w1)) 
}

###Our proposed method for estimating the unknown parameters, AUC, and ROC curve at given ss = (1:99)/100, along with their standard errors (SEs)###

ourmet<- function(y, R, x, v,ss)
{  
  
  
  para=estpara(y, R, x, v)
  mu=as.numeric(para$mu)
  phi=as.numeric(para$phi)
  separa=sqrt(as.numeric(diag(varpara(y,R,x,v,mu,phi))))
  k1=length(mu)
  outmu=cbind(mu,separa[1:k1])
  colnames(outmu)=c("Est","SE")
  
  outphi=cbind(phi,separa[c((k1+2):length(separa),k1+1)])
  colnames(outphi)=c("Est","SE")
  
  
  outauc=estauc(y,R,x,v,mu,phi)
  outauc=as.numeric(outauc)
  outauc=matrix(outauc,byrow=T,nrow=1)
  colnames(outauc)=c("Est","SE")
  rownames(outauc)=c("AUC")
  
  outroc=c()
  nss=length(ss)
  for(i in 1:nss)
  {  
  outroc1=estroc(ss[i],y,R,x,v,mu,phi)
  outroc1=as.numeric(outroc1)
  outroc=rbind(outroc,c(ss[i],outroc1))
  }
  
  
  colnames(outroc)=c("s","Est","SE")
  rownames(outroc)=1:nss
  
  
  
  
  list(outmu=outmu,outphi=outphi,outauc=outauc,outroc=outroc)  
}

##Proposed method for estimating the  AUC and ROC curve at given ss, along with their bootstrap SEs##

ourboot=function(y, R, x, v,ss,B=200)
{
 
set.seed(123456)
n=length(y)
out=ourest(y,R,x,v,ss)
out=as.numeric(out)
output=c()
for(i in 1:B)
{
bind=sample(1:n,n,replace=T)
by=y[bind]
bR=R[bind]
bx=x[bind]
bv=v[bind,]
bout=ourest(by,bR,bx,bv,ss)
bout=as.numeric(bout)
output=rbind(output,bout)
}  
bse=apply(output,2,sd)
bse=as.numeric(bse)

output=cbind(out,bse)
outauc=as.numeric(output[1,])
outauc=matrix(outauc,byrow=T,nrow=1)
colnames(outauc)=c("Est","SE")
rownames(outauc)=c("AUC")

nss=length(ss)
outroc=cbind(ss,output[-1,])
colnames(outroc)=c("s","Est","SE")
rownames(outroc)=1:nss

list(outauc=outauc,outroc=outroc)  

} 

##Goodness of fit test statistic for the verification model##
gofstat=function(y,R,x,v)
{  
  n=length(y)
  para=estpara(y, R, x, v)
  mu=as.numeric(para$mu)
  phi=as.numeric(para$phi)
  psi=phi[1:( length(phi)-1)]
  beta=phi[length(phi)]	
  xmat=cbind(rep(1,n),x, v)
  
  p1xv=1/(1+exp(as.numeric(xmat%*%mu)))
  cxv=log(1-p1xv+p1xv*exp(beta))
  txv=as.numeric(xmat%*%psi)+cxv
  pixv=1-plogis(txv)
  stat=sum((R-pixv)^2)-sum(pixv*(1-pixv))
  list(stat=stat,mu=mu,pixv=pixv)
}

##Goodness of fit test statistic and p-value for the verification model##
bootgof=function(y, R, x, v,B=200)
{
  

  n=length(x)
  out=gofstat(y,R,x,v)
  stat=out$stat
  mu=out$mu
  piest=out$pixv
  
  set.seed(123456)
  bootstat=rep(0,B)
  
  for(i in 1:B)
  {
    ind=sample(1:n,n,replace=T)
    newx=x[ind]
    newv=v[ind,]
    newR=rbinom(n,1,prob=piest[ind])
    newdismat1=cbind(rep(1,n),newx,newv)
    p1xv=1-plogis( as.numeric(newdismat1%*%mu) )
    newy=rbinom(n,1,p1xv)
    
    
    newtest=gofstat(newy,newR,newx,newv)$stat
    bootstat[i]=newtest
  }
  
  bse=sqrt(var(bootstat))
  tstat=stat/bse
  
  pvalue=2*(1-pnorm(abs(tstat)))
  output=c(stat,bse,tstat,pvalue)
  output=matrix(output,nrow=1)
  colnames(output)=c("T2","BSE","z-stat","p-value")
 output
  
}
