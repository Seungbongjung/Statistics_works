lda_cv=function(true_status,samp,partition){
  
  n=length(true_status); p=ncol(samp)  
  error_samp=numeric(length=partition);error_gl=numeric(length=partition)
  error_mpp=numeric(length=partition);error_map=numeric(length=partition)
  
  step.size=100
  tol=0.002
  P=matrix(1,p,p)
  diag(P)=0
  lam=0.06
  
  case=which(true_status==1); control=which(true_status==0)
  
  for(i in 1:partition){
    
    case_train=sort(sample(case,size=71,replace=FALSE))
    case_test=sort(setdiff(case,case_train))
    
    control_train=sort(sample(control,size=119,replace=FALSE))
    control_test=sort(setdiff(control,control_train))
    
    train=sort(c(case_train,control_train))
    test=sort(c(case_test,control_test))
    
    status_train=true_status[train]
    status_test=true_status[test]
    
    samp_train=samp[train,]
    samp_test=samp[test,]
    
    s1_train=samp[case_train,]; s0_train=samp[control_train,]
    
    mu1=numeric(length=p)
    for(j in 1:71){
      mu1=mu1+as.numeric(s1_train[j,])
    }
    mu1=mu1/71
    
    mu0=numeric(length=p)
    for(l in 1:119){
      mu0=mu0+as.numeric(s0_train[l,])
    }
    mu0=mu0/119
    
    est_cov_samp=1/190*t(samp_train)%*%(diag(190)-1/190*rep(1,190)%*%t(rep(1,190)))%*%samp_train+0.1^3*diag(p)
    est_prec_samp=solve(est_cov_samp)
    
    samp_cov=1/190*t(samp_train)%*%(diag(190)-1/190*rep(1,190)%*%t(rep(1,190)))%*%samp_train+0.001*diag(p)
    est_cov_gl=spcov(Sigma=samp_cov,S=samp_cov,lambda=lam*P,step.size=step.size,n.inner.steps=200,
                  thr.inner=0,tol.outer=tol,trace=1)$Sigma
    est_prec_gl=solve(est_cov_gl)
    
    mcmc_samp=mcmc(190,samp_cov,1,190*0.05,0.1,12000,3000,0.65)
    mcmc_mpp=mpp(mcmc_samp)
    mcmc_map=mmp(mcmc_samp)
    
    est_cov_mpp=bcd_covariance2(samp_cov,190,1,190*0.05,10000,0.1^3,mcmc_mpp)
    est_prec_mpp=solve(est_cov_mpp)
    
    est_cov_map=bcd_covariance2(samp_cov,190,1,190*0.05,10000,0.1^3,mcmc_map)
    est_prec_map=solve(est_cov_map)
    
    p1=71/190
    p0=119/190
    
    error_temp_samp=0
    error_temp_gl=0
    error_temp_proposed_mpp=0
    error_temp_proposed_map=0
    
    for(ind in 1:378){
      est_samp=0; est_gl=0l est_mpp=0; est_map=0
      
      ob=samp_test[ind,]
      rule1_samp=as.numeric(t(ob)%*%est_prec_samp%*%mu1-1/2*t(mu1)%*%est_prec_samp%*%mu1+log(p1))
      rule0_samp=as.numeric(t(ob)%*%est_prec_samp%*%mu0-1/2*t(mu0)%*%est_prec_samp%*%mu0+log(p0))
      
      rule1_gl=as.numeric(t(ob)%*%est_prec_gl%*%mu1-1/2*t(mu1)%*%est_prec_gl%*%mu1+log(p1))
      rule0_gl=as.numeric(t(ob)%*%est_prec_gl%*%mu0-1/2*t(mu0)%*%est_prec_gl%*%mu0+log(p0))
      
      rule1_mpp=as.numeric(t(ob)%*%est_prec_mpp%*%mu1-1/2*t(mu1)%*%est_prec_mpp%*%mu1+log(p1))
      rule0_mpp=as.numeric(t(ob)%*%est_prec_mpp%*%mu0-1/2*t(mu0)%*%est_prec_mpp%*%mu0+log(p0))
      
      rule1_map=as.numeric(t(ob)%*%est_prec_map%*%mu1-1/2*t(mu1)%*%est_prec_map%*%mu1+log(p1))
      rule0_map=as.numeric(t(ob)%*%est_prec_map%*%mu0-1/2*t(mu0)%*%est_prec_map%*%mu0+log(p0))
      
      if(rule1_samp>rule0_samp){
        est_samp=1
      }
      if(rule1_gl>rule0_gl){
        est_gl=1
      }
      if(rule1_mpp>rule0_mpp){
        est_mpp=1
      }
      if(rule1_map>rule0_map){
        est_map=1
      }
      
      diff_samp=est_samp-status_test[ind]
      error_temp_samp=error_temp_samp+abs(diff)
    
      diff_gl=est_gl-status_test[ind]
      error_temp_gl=error_temp_gl+abs(diff)
      
      
      diff_mpp=est_mpp-status_test[ind]
      error_mpp_samp=error_temp_mpp+abs(diff)

      diff_map=est_map-status_test[ind]
      error_temp_map=error_temp_map+abs(diff)
    }
    error_samp[i]=error_temp_samp/378
    error_gl[i]=error_temp_gl/378
    error_mpp[i]=error_temp_mpp/378
    error_map[i]=error_temp_map/378
  }
  error=list(error_samp,error_gl,error_mpp,error_map)
  return(error)
}