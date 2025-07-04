rm(list=ls())
source("gendata.R") 
source("supp.R")
source("Our.R")
source("Others.R")
rs=read.table("seed.txt")[,1]

###Generate data from the first random seed under Scenario 5 with n=1000###
set.seed(rs[1])
n=1000
dat=simudata5(n)
y=dat$y
R=dat$R
v=dat$v
x=dat$x

###Our proposed method for estimating the unknown parameters, AUC, and ROC curve at given ss = (1:99)/100, along with their standard errors (SEs)###
### ###
ss=(1:99)/100
out=ourmet(y,R,x,v,ss)

##mu estimates and SEs##
out$outmu

##phi estimates and SES##
out$outphi

##AUC estimate and SE##
out$outauc

##ROC(s) estimates and SEs##
out$outroc

##Proposed method for estimating the  AUC and ROC curve at given ss = (1:99)/100, along with their bootstrap SEs##
##Number of bootstrap sample is B=200##
##It may take around 2 to 5 minutes to run the code##
ss=(1:99)/100
outb=ourboot(y,R,x,v,ss,B=200)

outb$outauc
outb$outroc

###Goodness of fit test for the disease model###
###For this example: p-value=0.6930; there is no evidence to reject the logistic disease model###
library(rms)
yp=y[R==1]
xp=x[R==1]
vp=v[R==1,]
outp=lrm(yp~xp+vp,x=T,y=T)
round(resid(outp,"gof"),4)

##Goodness of fit test for the verification model##
##It may take around 2 to 5 minutes to run the code##
######For this example: p-value=0.962; there is no evidence to reject the logistic verification model###

bootgof(y,R,x,v,B=200)

###IG method: point estimates of AUC and ROC at ss=(1:99)/100###
ss=(1:99)/100
igmet(y, R, x, v,ss)

###VER method: point estimates of AUC and ROC at ss=(1:99)/100###
ss=(1:99)/100
vermet(y, R, x, ss)

###Full method: point estimates of AUC and ROC at ss=(1:99)/100###
ss=(1:99)/100
fulmet(y,  x,ss)

###IPW method: point estimate of AUC###

source("supp.R")
source("IPW.R")
ipwest(y,R,x,v)
  
###IPW method: point estimate of AUC and SE of the AUC estimate###
source("supp.R")
source("IPW.R")
library("Rcpp")
sourceCpp("ipw.cpp")

ipwmet(y,R,x,v)
