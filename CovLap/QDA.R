qda_samp=function(true_status,samp){
  
  n=length(true_status); p=ncol(samp)  
  error=0
  
  where=NULL

  for(ind in 1:n){

    samp_temp=samp[-ind,]
    true_temp=true_status[-ind]
    n_temp=n-1
    
    u1=which(true_temp==1); u0=which(true_temp==0)
    l1=length(u1); l0=length(u0)
    p1=l1/n_temp; p0=1-p1
    
    s1=samp_temp[u1,];  s0=samp_temp[u0,]
    
    mu1=numeric(length=p)
    for(j in 1:length(u1)){
      mu1=mu1+as.numeric(s1[j,])
    }
    mu1=mu1/l1
    
    mu0=numeric(length=p)
    for(l in 1:length(u0)){
      mu0=mu0+as.numeric(s0[l,])
    }
    mu0=mu0/l0
  
    est_cov0=1/l0*t(s0)%*%(diag(l0)-1/l0*rep(1,l0)%*%t(rep(1,l0)))%*%s0+0.001*diag(p)
    est_prec0=solve(est_cov0)
    
    est_cov1=1/l1*t(s1)%*%(diag(l1)-1/l1*rep(1,l1)%*%t(rep(1,l1)))%*%s1+0.001*diag(p)
    est_prec1=solve(est_cov1)
  
    est=0
    ob=samp[ind,]
    rule1=as.numeric(-1/2*log(det(est_cov1))-1/2*t(ob-mu1)%*%est_prec1%*%(ob-mu1)+log(p1))
    rule0=as.numeric(-1/2*log(det(est_cov0))-1/2*t(ob-mu0)%*%est_prec0%*%(ob-mu0)+log(p0))
    
    if(rule1>rule0){
      est=1
    }
    diff=est-true_status[ind]
    if(diff!=0){
      where=c(where,ind)
    }
    error=error+abs(diff)
    
  }
  error=error/n
  return(list(error,where))
}

qda_gl=function(true_status,samp){
  
  p=ncol(samp)  
  
  step.size=100
  tol=0.001
  P=matrix(1,p,p)
  diag(P)=0
  lam=0.06
  
  n=length(true_status)
  error=0
  
  where=NULL
  
  for(ind in 1:n){
    
    print(paste(ind,"th sample",sep=""))
    
    samp_temp=samp[-ind,]
    true_temp=true_status[-ind]
    n_temp=n-1
    
    u1=which(true_temp==1); u0=which(true_temp==0)
    l1=length(u1); l0=length(u0)
    
    p1=l1/n_temp; p0=1-p1
    
    s1=samp_temp[u1,];  s0=samp_temp[u0,]
    
    mu1=numeric(length=p)
    for(j in 1:l1){
      mu1=mu1+as.numeric(s1[j,])
    }
    mu1=mu1/l1
    
    mu0=numeric(length=p)
    for(l in 1:l0){
      mu0=mu0+as.numeric(s0[l,])
    }
    mu0=mu0/l0
    
    samp_cov0=1/l0*t(s0)%*%(diag(l0)-1/l0*rep(1,l0)%*%t(rep(1,l0)))%*%s0+0.001*diag(p)
    est_cov0=spcov(Sigma=samp_cov0,S=samp_cov0,lambda=lam*P,step.size=step.size,n.inner.steps=200,
          thr.inner=0,tol.outer=tol,trace=1)$Sigma
    est_prec0=solve(est_cov0)
    
    samp_cov1=1/l1*t(s1)%*%(diag(l1)-1/l1*rep(1,l1)%*%t(rep(1,l1)))%*%s1+0.001*diag(p)
    est_cov1=spcov(Sigma=samp_cov1,S=samp_cov1,lambda=lam*P,step.size=step.size,n.inner.steps=200,
                   thr.inner=0,tol.outer=tol,trace=1)$Sigma
    est_prec1=solve(est_cov1)
    
    est=0
    ob=samp[ind,]
    rule1=as.numeric(-1/2*log(det(est_prec1))-1/2*t(ob-mu1)%*%est_prec1%*%(ob-mu1)+log(p1))
    rule0=as.numeric(-1/2*log(det(est_prec0))-1/2*t(ob-mu0)%*%est_prec0%*%(ob-mu0)+log(p0))
    
    if(rule1>rule0){
      est=1
    }
    diff=est-true_status[ind]
    if(diff!=0){
      where=c(where,ind)
    }
    error=error+abs(diff)
  }
  error=error/n
  return(list(error,where))
}

qda_proposed=function(true_status,samp){
  
  n=length(true_status); p=ncol(samp)  
  error_mpp=0; error_map=0
  
  for(ind in 1:40){
    
    n_temp=n-1
    u1=which(true_temp==1); u0=which(true_temp==0)
    l1=length(u1); l0=length(u0)
    p1=l1/n_temp; p0=1-p1
    
    s1=samp_temp[u1,];  s0=samp_temp[u0,]
    
    mu1=numeric(length=p)
    for(j in 1:length(u1)){
      mu1=mu1+as.numeric(s1[j,])
    }
    mu1=mu1/l1
    
    mu0=numeric(length=p)
    for(l in 1:length(u0)){
      mu0=mu0+as.numeric(s0[l,])
    }
    mu0=mu0/l0
    
    samp_cov0=1/l0*t(s0)%*%(diag(l0)-1/l0*rep(1,l0)%*%t(rep(1,l0)))%*%s0+0.001*diag(p)
    mcmc_samp0=mcmc(l0,samp_cov0,1,l0*0.06,0.1^3,5000,5000,4/(p-1))
    mcmc_mpp0=mpp(mcmc_samp0)
    mcmc_map0=mmp(mcmc_samp0)
    
    est_cov_mpp0=bcd_covariance2(samp_cov0,l0,1,l0*0.06,10000,0.1^3,mcmc_mpp0)
    est_prec_mpp0=solve(est_cov_mpp0)
    
    est_cov_map0=bcd_covariance2(samp_cov0,l0,1,l0*0.06,10000,0.1^3,mcmc_map0)
    est_prec_map0=solve(est_cov_map0)
    
    samp_cov1=1/l1*t(s1)%*%(diag(l1)-1/l1*rep(1,l1)%*%t(rep(1,l1)))%*%s1+0.001*diag(p)
    mcmc_samp1=mcmc(l1,samp_cov1,1,l1*0.06,0.1^3,5000,5000,4/(p-1))
    mcmc_mpp1=mpp(mcmc_samp1)
    mcmc_map1=mmp(mcmc_samp1)
    
    est_cov_mpp1=bcd_covariance2(samp_cov1,l1,1,l1*0.06,10000,0.1^3,mcmc_mpp1)
    est_prec_mpp1=solve(est_cov_mpp1)
    
    est_cov_map1=bcd_covariance2(samp_cov1,l1,1,l1*0.06,10000,0.1^3,mcmc_map1)
    est_prec_map1=solve(est_cov_map1)
    
    est_mpp=0; est_map=0
    ob=samp[ind,]
    rule1_mpp=as.numeric(-1/2*log(det(est_prec_mpp1))-1/2*t(ob-mu1)%*%est_prec_mpp1%*%(ob-mu1)+log(p1))
    rule0_mpp=as.numeric(-1/2*log(det(est_prec_mpp0))-1/2*t(ob-mu0)%*%est_prec_mpp0%*%(ob-mu0)+log(p0))
  
    if(rule1_mpp>rule0_mpp){
      est_mpp=1
    }
    
    rule1_map=as.numeric(-1/2*log(det(est_prec_map1))-1/2*t(ob-mu1)%*%est_prec_map1%*%(ob-mu1)+log(p1))
    rule0_map=as.numeric(-1/2*log(det(est_prec_map0))-1/2*t(ob-mu0)%*%est_prec_map0%*%(ob-mu0)+log(p0))
    if(rule1_map>rule0_map){
      est_map=1
    }
    diff_mpp=est_mpp-true_status[ind]
    error_mpp=error_mpp+abs(diff_mpp)
    
    diff_map=est_map-true_status[ind]
    error_map=error_map+abs(diff_map)
    print(paste("mpp error: ",error_mpp," map error: ",error_map,sep=""))
  }
  return(c(error_mpp,error_map))
}