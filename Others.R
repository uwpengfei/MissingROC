############################ IG method##### 
igmet <- function(y, R, x, v,ss)
{
  n=length(y)
  xmat=cbind(x,v)
  xmat=as.matrix(xmat)
  dismat=xmat ##Covariates for disease model
  
  mu0=-glm(y[R>0]~dismat[R>0,],family="binomial")$coef
  mu=as.numeric(mu0)
  mub=c(mu,0)
  gmat=cbind(rep(1,n),xmat,(R-1))
  gmat=as.matrix(gmat)
  gxvr=1/( 1 + exp(as.numeric(gmat%*%mub)) )
  w1=(gxvr/sum(gxvr))
  w0=((1-gxvr)/sum(1-gxvr))
  
  outauc=newauc(x,w0,w1)
  outroc=cbind(ss,newroc(ss,x,w0,w1))
  colnames(outroc)=c("s","Est")
  rownames(outroc)=1:length(ss)
  
list(outauc=outauc,outroc=outroc)
}

############### VER method#####

vermet=function(y,R,x,ss)
{
  y=y[R==1]
  x=x[R==1]
  w0=1-y
  w1=y
  
  outauc=newauc(x,w0,w1)
  outroc=cbind(ss,newroc(ss,x,w0,w1))
  colnames(outroc)=c("s","Est")
  rownames(outroc)=1:length(ss)
  
  list(outauc=outauc,outroc=outroc)

}






############### Full method #####

fulmet=function(y,x,ss)
{
  x1=x[y==1]
  x0=x[y==0]	
  w0=1-y
  w1=y
  
  outauc=newauc(x,w0,w1)
  outroc=cbind(ss,newroc(ss,x,w0,w1))
  colnames(outroc)=c("s","Est")
  rownames(outroc)=1:length(ss)
  
  list(outauc=outauc,outroc=outroc) 

}
















