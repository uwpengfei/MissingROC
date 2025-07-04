##ewcdf: calcualte the cdf at distinct values##
ewcdf=function (x, w ) 
{
  w=w/sum(w)
  tt=sort(unique(x))
  obj=function(t)
  {
    sum(w[x<=t])    
  }  
  
  estF=as.numeric(sapply(tt,obj))
  estF
}

## a more efficient way to estimate auc
newauc=function(x,w0,w1)
{
  w0=w0/sum(w0)
  w1=w1/sum(w1)
  hatF0=ewcdf(x,w0)
  hatF1=ewcdf(x,w1)
  neww0=c(hatF0[1],diff(hatF0))
  neww1=c(hatF1[1],diff(hatF1))
  
  out1=sum(neww1*hatF0)
  out2=0.5*sum(neww0*neww1)
  out1-out2
}

###################### ROC curve estimation
disroc=function(s,hatF0,hatF1)
{
  n=length(hatF0)
  
  ind=min(c((1:n)[hatF0>=1-s],n))
  1-hatF1[ind]
}

newroc=function(ss,x,w0,w1)
{
estF0=ewcdf(x,w0)
estF1=ewcdf(x,w1)

estroc=as.numeric(sapply(ss,disroc,hatF0=estF0,hatF1=estF1))
estroc
}  






