
library("nleqslv")
library("MASS")

######################### IPW method#####
mu_est <- function(y, R, cova){
  yr <- ifelse(R == 1, y, -1) 
  # use glm, logistic regression directly
  ind = which(R == 1)
  cova.obs = cova[ind, ]
  yobs = yr[ind]
  
  result = glm(yobs ~ cova.obs , family="binomial") 
  mu <- -coef(result)  
  
  mu[is.na(mu)] = 0  ## avoid constant covariate
  return(mu)
}


## estimate parameters phi in verification model using nleqslv package ###
phi_EM_nleqslv <- function(init_phi, mu, y, R, cova1, cova2, em_run = 1000){
  yr <- ifelse(R == 1, y, -1) 
  n = length(R)
  
  k1 = ncol(cova1)
  k2 = ncol(cova2)
  
  old_phi = init_phi
  
  # objective score function
  sbarfun <- function(z){
    sfun <- function(ty){
      s = matrix(0, n, k2 + 2)
      temps = temp0 + z[k2 + 2] * ty
      
      tps = 1/(1 + exp(temps))
      tt = tps - R
      s[, 1] = tt
      s[, k2+2] = ty * tt
      
      for(i in 1:k2){
        s[, i + 1] =  cova2[, i] * tt
      }
      
      return(s)
    }
    temp0 = z[1]
    for(i in 1:k2) temp0 =  temp0 + z[1 + i] * cova2[, i]
    
    s1 = sfun(1)
    s0 = sfun(0)
    E0s = p0 * s1 + (1 - p0) * s0
    sbar = colSums(R * sfun(yr)) + colSums((1 - R) * E0s)
    return(sbar/n)
  }
  
  # jacobian for the objective function
  jac_fun <- function(z){
    DS2 = matrix(0, k2+2, k2+2)
    b0 = c(rep(0, k2 + 1), 1)
    pi_hat = rep(0, n)
    for(i in 1:n){
      ay = c(1, cova2[i, ], yr[i])
      a1 = c(1, cova2[i, ], 1)
      a0 = c(1, cova2[i, ], 0)
      pi1 = 1/(1 + exp(sum(z * a1)))      ## p(R=1|y=1)
      pi0 = 1/(1 + exp(sum(z * a0)))      ## p(R=1|y=0)
      
      pi_hat[i] = 1/(1 + exp(sum(z * ay)))
      ds1 = -pi1 * (1 - pi1) * outer(a1, a1)
      ds0 = -pi0 * (1 - pi0) * outer(a0, a0)
      DS2 = DS2 - R[i] * pi_hat[i] * (1 - pi_hat[i]) * outer(ay, ay) + 
        (1 - R[i]) * (p0[i] * (1 - p0[i]) * outer((pi1 - pi0) * a0 + (pi1 - R[i]) * b0, b0 ) +
                        p0[i] * (ds1 - ds0) + ds0 ) 
    }
    return(DS2)
  }
  
  temp = mu[1]
  for(i in 1:k1){
    temp = temp + mu[1+i] * cova1[, i]
  }
  p1 = 1/(1 + exp(temp))   ## p1=pr(y=1|R=1, x, v)
  
  # using EM
  run = 0
  repeat{
    run = run + 1
    
    p0 = p1 * exp(old_phi[k2+2])/(1 - p1 + p1 * exp(old_phi[k2+2]))
    
    result = nleqslv(old_phi, fn = sbarfun) 
    phi <- result$x
    if(max(abs(old_phi - phi)) < 0.0001 || run == em_run) break
    
    old_phi = phi
    
  }
  
  return(list('phi' = phi, 'fval' = sum(abs(result$fvec)), 'convN' = run))
}




###IPW method: variance estimate of the AUC estimate###
var_auc_iv <- function(pi_hat, mu, phi, yr, R, x, v, auc1){ 
  n = length(R)
  mu=as.numeric(mu)
  phi=as.numeric(phi)
  auc1=as.numeric(auc1)
  k=length(mu)
  
  cova=cbind(rep(1,n),x,v)
  cova2=cbind(x,v)
  
  p1 = 1/(1 + exp(cova%*%mu))    ## pr(y=1|R=1)
  beta=phi[k+1]
  p0 = p1 * exp(beta)/(1 - p1 + p1 * exp(beta))  ## pr(y=1|R=0) 
  
  ## calculate the first term of variance, using c
  F1 = cal_F1(pi_hat, auc1, yr, x, R)
  
  # calculate Gamma1 and Gamma2, using c
  
  gamma1 = cal_gamma1_v2(pi_hat, auc1, yr, R, cova[,-1])
  
  # calculate s2 and derivative of s2
  s2 = matrix(0, n, k+1)
  s1 = matrix(0, n, k)
  FIM2 = matrix(0, k+1, k+1)
  FIM1 = matrix(0, k, k)
  b0 = c(rep(0, k), 1)
  ds2_mu = matrix(0, k, k+1)
  for(i in 1:n){
    ay = c(1, cova2[i, ], yr[i])
    a1 = c(1, cova2[i, ], 1)
    a0 = c(1, cova2[i, ], 0)
    pi1 = 1/(1 + exp(sum(phi * a1)))      ## p(R=1|y=1)
    pi0 = 1/(1 + exp(sum(phi * a0)))      ## p(R=1|y=0)
    s2[i, ] = R[i] * (pi_hat[i] - R[i]) * ay + (1 - R[i]) * 
      (p0[i] * (pi1 - R[i]) * a1 + (1 - p0[i]) * (pi0 - R[i]) * a0)
    
    
    ds2_1 = -pi1 * (1 - pi1) * outer(a1, a1)
    ds2_0 = -pi0 * (1 - pi0) * outer(a0, a0)
    FIM2 = FIM2 - R[i] * pi_hat[i] * (1 - pi_hat[i]) * outer(ay, ay) + 
      (1 - R[i]) * (p0[i] * (1 - p0[i]) * outer((pi1 - pi0) * a0 + (pi1 - R[i]) * b0, b0 ) +
                      p0[i] * (ds2_1 - ds2_0) + ds2_0 ) 
    
    ## add variance term comes from estimation mu (added since 10/23/2016)
    ax = cova[i, ]
    s1[i, ] = R[i] * (p1[i] - yr[i]) * ax 
    FIM1 = FIM1 - R[i] * p1[i] * (1 - p1[i]) * outer(ax, ax)
    ds2_mu = ds2_mu - (1 - R[i]) * p1[i] * (1 - p1[i]) * exp(beta)/(1-p1[i]+p1[i]*exp(beta))^2 *
      (pi1 * outer(ax, a1) - pi0 * outer(ax, a0)) 
  }
  ds2_mu = ds2_mu/n
  FIM1 = FIM1/n
  FIM2 = FIM2/n
  gamma1 = gamma1/n/(n-1)
  
  
  vq = var(F1 - (s2 - s1 %*% ginv(FIM1) %*% ds2_mu) %*% ginv(FIM2) %*% gamma1)
  
  dp1 = sum(R * yr / pi_hat)/n    ## estimation of pr(y=1)
  dp0 = sum(R * (1 - yr)/ pi_hat)/n  ## another way of estimating pr(y=0)
  #dp0 = 1- dp1
  
  vA = vq/(dp1 * dp0) ^2
  
  
  return(vA/n)
}



###Point estimates of unknown parameters and AUC###
ipwest=function(y,R,x,v)
  {
  n=length(y)
  yr <- ifelse(R == 1, y, -1) 
  cova = cbind(x, v)
  cova1 = as.matrix(cova)
  cova2 = as.matrix(cova)	
  mu = mu_est(yr, R, cova1)
  k2=ncol(cova2)+1
  k1=ncol(cova1)+1
  
  calauc=function(para)
  {
    phi=para[(k1+1):(k1+k2+1)]
    
    pimat=cbind(rep(1,n),cova2,y)
    pixvy=1/( 1 + exp(as.numeric(pimat%*%phi)) )
    w1=y*R/pixvy
    w0=(1-y)*R/pixvy
    
    prauc=newauc(x,w0,w1)
    prauc
  }  
  
  init_phi = c(1,rep(0, k2) )
  
  temp_res = phi_EM_nleqslv(init_phi, mu, yr, R, cova1, cova2)
  
  while(temp_res$fval>1e-06) 
  {
    init_phi = c(1,rep(0, k2 ) )+runif(k2+1,-0.5,0.5)
    temp_res = phi_EM_nleqslv(init_phi, mu, yr, R, cova1, cova2)
  }
  
  mu=mu
  phi=temp_res$phi
  
  par= c(mu,temp_res$phi)
  
  
  phi=temp_res$phi
  AUC=calauc(par)
  
  
  list(phi=phi,mu=mu,AUC=AUC)
}

###Point estimates of unknown parameters and AUC, along with the SE of the AUC estimate###
ipwmet=function(y,R,x,v)
  {
  out=ipwest(y,R,x,v)
  mu=out$mu
  phi=out$phi
  AUC=out$AUC
  
  xxmat=cbind(rep(1,n),x,v,y)
  xxmat=as.matrix(xxmat)
  pihat=1/(1+exp(xxmat%*%as.numeric(phi)))
  variance=var_auc_iv(pi_hat=pihat, mu=mu, phi=phi, yr=ifelse(R==1,y,-1),R=R,x=x,v=v,auc1 = AUC)
  outest=cbind(AUC,as.numeric(sqrt(variance)) )
  colnames(outest)=c("Est","SE")
  
  list(phi=phi,mu=mu,outest=outest)
}






