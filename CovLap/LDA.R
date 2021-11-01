library(spcov)

lda_samp=function(true_status,samp){
  
  n=length(true_status); p=ncol(samp)  
  error=0
  
  where=NULL

  for(ind in 1:n){

    samp_temp=samp[-ind,]
    true_temp=true_status[-ind]
    n_temp=n-1
    
    u1=which(true_temp==1); u0=which(true_temp==0)
    p1=length(u1)/n_temp; p0=1-p1
    
    s1=samp_temp[u1,];  s0=samp_temp[u0,]
    
    mu1=numeric(length=p)
    for(j in 1:length(u1)){
      mu1=mu1+as.numeric(s1[j,])
    }
    mu1=mu1/length(u1)
    
    mu0=numeric(length=p)
    for(l in 1:length(u0)){
      mu0=mu0+as.numeric(s0[l,])
    }
    mu0=mu0/length(u0)
    
    est_cov=1/n_temp*t(samp_temp)%*%(diag(n_temp)-1/n_temp*rep(1,n_temp)%*%t(rep(1,n_temp)))%*%samp_temp
    est_prec=solve(est_cov)
  
    est=0
    ob=samp[ind,]
    rule1=as.numeric(t(ob)%*%est_prec%*%mu1-1/2*t(mu1)%*%est_prec%*%mu1+log(p1))
    rule0=as.numeric(t(ob)%*%est_prec%*%mu0-1/2*t(mu0)%*%est_prec%*%mu0+log(p0))
    
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

lda_gl=function(true_status,samp){
  
  step.size=100
  tol=0.001
  P=matrix(1,p,p)
  diag(P)=0
  lam=0.06
  
  n=length(true_status); p=ncol(samp)  
  error=0
  
  for(ind in 1:n){
    
    samp_temp=samp[-ind,]
    true_temp=true_status[-ind]
    n_temp=n-1
    
    u1=which(true_temp==1); u0=which(true_temp==0)
    p1=length(u1)/n_temp; p0=1-p1
    
    s1=samp_temp[u1,];  s0=samp_temp[u0,]
    
    mu1=numeric(length=p)
    for(j in 1:length(u1)){
      mu1=mu1+as.numeric(s1[j,])
    }
    mu1=mu1/length(u1)
    
    mu0=numeric(length=p)
    for(l in 1:length(u0)){
      mu0=mu0+as.numeric(s0[l,])
    }
    mu0=mu0/length(u0)
    
    samp_cov=1/n_temp*t(samp_temp)%*%(diag(n_temp)-1/n_temp*rep(1,n_temp)%*%t(rep(1,n_temp)))%*%samp_temp
    samp_cov=samp_cov+0.001*diag(p)
    
    est_cov=spcov(Sigma=samp_cov,S=samp_cov,lambda=lam*P,step.size=step.size,n.inner.steps=200,
                  thr.inner=0,tol.outer=tol,trace=1)$Sigma
    est_prec=solve(est_cov)
    
    est=0
    ob=samp[ind,]
    rule1=as.numeric(t(ob)%*%est_prec%*%mu1-1/2*t(mu1)%*%est_prec%*%mu1+log(p1))
    rule0=as.numeric(t(ob)%*%est_prec%*%mu0-1/2*t(mu0)%*%est_prec%*%mu0+log(p0))
    
    if(rule1>rule0){
      est=1
    }
    diff=est-true_status[ind]
    error=error+abs(diff)
    
  }
  error=error/n
  return(error)
}

lda_proposed=function(true_status,samp){
  
  n=length(true_status); p=ncol(samp)  
  error_mpp=0; error_map=0
  
  for(ind in 1:40){
    
    print(paste(c("replication: ",ind),sep=""))
    samp_temp=samp[-ind,]
    true_temp=true_status[-ind]
    n_temp=n-1
    
    u1=which(true_temp==1); u0=which(true_temp==0)
    p1=length(u1)/n_temp; p0=1-p1
    
    s1=samp_temp[u1,];  s0=samp_temp[u0,]
    
    mu1=numeric(length=p)
    for(j in 1:length(u1)){
      mu1=mu1+as.numeric(s1[j,])
    }
    mu1=mu1/length(u1)
    
    mu0=numeric(length=p)
    for(l in 1:length(u0)){
      mu0=mu0+as.numeric(s0[l,])
    }
    mu0=mu0/length(u0)
    
    samp_cov=1/n_temp*t(samp_temp)%*%(diag(n_temp)-1/n_temp*rep(1,n_temp)%*%t(rep(1,n_temp)))%*%samp_temp
    samp_cov=samp_cov+0.001*diag(p)
    sum(edge_ind(samp_cov,0.1^3))
    
    mcmc_samp=mcmc(n_temp,samp_cov,1,n_temp*0.06,0.1^3,5000,5000,4/21)
    mcmc_mpp=mpp(mcmc_samp)
    mcmc_map=mmp(mcmc_samp)
    
    est_cov_mpp=bcd_covariance2(samp_cov,n_temp,1,n_temp*0.06,10000,0.1^3,mcmc_mpp)
    est_prec_mpp=solve(est_cov_mpp)
    
    est_cov_map=bcd_covariance2(samp_cov,n_temp,1,n_temp*0.06,10000,0.1^3,mcmc_map)
    est_prec_map=solve(est_cov_map)
    
    est_mpp=0; est_map=0
    ob=samp[ind,]
    rule1_mpp=as.numeric(t(ob)%*%est_prec_mpp%*%mu1-1/2*t(mu1)%*%est_prec_mpp%*%mu1+log(p1))
    rule0_mpp=as.numeric(t(ob)%*%est_prec_mpp%*%mu0-1/2*t(mu0)%*%est_prec_mpp%*%mu0+log(p0))
  
    if(rule1_mpp>rule0_mpp){
      est_mpp=1
    }
    
    rule1_map=as.numeric(t(ob)%*%est_prec_map%*%mu1-1/2*t(mu1)%*%est_prec_map%*%mu1+log(p1))
    rule0_map=as.numeric(t(ob)%*%est_prec_map%*%mu0-1/2*t(mu0)%*%est_prec_map%*%mu0+log(p0))
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