################ code for generating data  from Scenario 1########################
##Input: n, sample size##
##Output: y, R, x, v##

simudata1=function(n)
{

  
  mu=c(1.7,-1.5,-1.5,-2.5)
  phi=c(1.3,-1.2,1,-1.5,0)
  
  v1=rnorm(n,0,1)
  v2=rbinom(n,1,0.5)
  x=runif(n,-1,1)
  xmat=cbind(rep(1,n),v1,v2,x)
  
  beta=phi[length(phi)]
  phi0=phi[-(length(phi))]
  
  disprob1=1-plogis(as.numeric(xmat%*%mu) )
  cxv=log( 1-disprob1+disprob1*exp(beta) )
  txv=as.numeric(xmat%*%phi0+cxv)
  
  misprob=1-plogis(txv)
  R=rbinom(n,1,misprob)
  
  disprob=1-plogis(as.numeric(xmat%*%mu+(R-1)*beta ) )
  y=rbinom(n,1,disprob)
  
  list(y=y,R=R,x=x,v=cbind(v1,v2)) 
}  


################ code for generating data  from Scenario 2########################
##Input: n, sample size##
##Output: y, R, x, v##

simudata2=function(n)
{
  ##X: biomarker
  ##V1 and V2: two covariates
  mu=c(1.7,-1.5,-1.5,-2.5)
  phi=c(1.3,-1.2,1,-1.5,-2)
  
  v1=rnorm(n,0,1)
  v2=rbinom(n,1,0.5)
  x=runif(n,-1,1)
  xmat=cbind(rep(1,n),v1,v2,x)
  
  beta=phi[length(phi)]
  phi0=phi[-(length(phi))]
  
  disprob1=1-plogis(as.numeric(xmat%*%mu) )
  cxv=log( 1-disprob1+disprob1*exp(beta) )
  txv=as.numeric(xmat%*%phi0+cxv)
  
  misprob=1-plogis(txv)
  R=rbinom(n,1,misprob)
  
  disprob=1-plogis(as.numeric(xmat%*%mu+(R-1)*beta ) )
  y=rbinom(n,1,disprob)
  
  list(y=y,R=R,x=x,v=cbind(v1,v2))  
}  

################ code for generating data  from Scenario 3########################
##Input: n, sample size##
##Output: y, R, x, v##

simudata3=function(n)
{
  ##X: biomarker
  ##V1 and V2: two covariates
  mu=c(1.7,-1.5,0.5,-1.5,-2.5)
  phi=c(1.3,-1.2,0.5,1,-1.5,-2)
  
  v1=rnorm(n,0,1)
  v2=rbinom(n,1,0.5)
  x=runif(n,-1,1)
  xmat=cbind(rep(1,n),v1,v1^2,v2,x)
  
  beta=phi[length(phi)]
  phi0=phi[-(length(phi))]
  
  disprob1=1-plogis(as.numeric(xmat%*%mu) )
  cxv=log( 1-disprob1+disprob1*exp(beta) )
  txv=as.numeric(xmat%*%phi0+cxv)
  
  misprob=1-plogis(txv)
  R=rbinom(n,1,misprob)
  
  disprob=1-plogis(as.numeric(xmat%*%mu+(R-1)*beta ) )
  y=rbinom(n,1,disprob)
  
  list(y=y,R=R,x=x,v=cbind(v1,v2))   
}  


################ code for generating data  from Scenario 4########################
##Input: n, sample size##
##Output: y, R, x, v##

simudata4=function(n)
{
  ##X: biomarker
  ##V1 and V2: two covariates
  mu=c(1.7,-1.5,-1.5,-2.5)
  phi=-c(1.3,-1.5,1,-1.5,0)
  
  v1=rnorm(n,0,1)
  v2=rbinom(n,1,0.5)
  x=runif(n,-1,1)
  xmat=cbind(rep(1,n),v1,v2,x)
  
  beta=phi[length(phi)]
  phi0=phi[-(length(phi))]
  
  disprob1=1-plogis(as.numeric(xmat%*%mu) )
  cxv=log( 1-disprob1+disprob1*exp(beta) )
  txv=as.numeric(xmat%*%phi0+cxv)
  
  misprob=1-plogis(txv)
  R=rbinom(n,1,misprob)
  
  disprob=1-plogis(as.numeric(xmat%*%mu+(R-1)*beta ) )
  y=rbinom(n,1,disprob)
  
  list(y=y,R=R,x=x,v=cbind(v1,v2)) 
}  


################ code for generating data  from Scenario 5########################
##Input: n, sample size##
##Output: y, R, x, v##

simudata5=function(n)
{
  ##X: biomarker
  ##V1 and V2: two covariates
  mu=c(1.7,-1.5,-1.5,-3)
  phi=-c(0.5,-1.5,1,-1.5,3)
  
  v1=rnorm(n,0,1)
  v2=rbinom(n,1,0.5)
  x=runif(n,-1,1)
  xmat=cbind(rep(1,n),v1,v2,x)
  
  beta=phi[length(phi)]
  phi0=phi[-(length(phi))]
  
  disprob1=1-plogis(as.numeric(xmat%*%mu) )
  cxv=log( 1-disprob1+disprob1*exp(beta) )
  txv=as.numeric(xmat%*%phi0+cxv)
  
  misprob=1-plogis(txv)
  R=rbinom(n,1,misprob)
  
  disprob=1-plogis(as.numeric(xmat%*%mu+(R-1)*beta ) )
  y=rbinom(n,1,disprob)
  
  list(y=y,R=R,x=x,v=cbind(v1,v2))  
}  


################ code for generating data  from Scenario 6########################
##Input: n, sample size##
##Output: y, R, x, v##


simudata6=function(n)
{
  ##X: biomarker
  ##V1 and V2: two covariates
  mu=c(1.7,-1.5,0.5,-1.5,-3)
  phi=-c(0.5,-1.5,0.5,1,-1.5,3)
  
  v1=rnorm(n,0,1)
  v2=rbinom(n,1,0.5)
  x=runif(n,-1,1)
  xmat=cbind(rep(1,n),v1,v1^2,v2,x)
  
  beta=phi[length(phi)]
  phi0=phi[-(length(phi))]
  
  disprob1=1-plogis(as.numeric(xmat%*%mu) )
  cxv=log( 1-disprob1+disprob1*exp(beta) )
  txv=as.numeric(xmat%*%phi0+cxv)
  
  misprob=1-plogis(txv)
  R=rbinom(n,1,misprob)
  
  disprob=1-plogis(as.numeric(xmat%*%mu+(R-1)*beta ) )
  y=rbinom(n,1,disprob)
  
  list(y=y,R=R,x=x,v=cbind(v1,v2))   
}  



