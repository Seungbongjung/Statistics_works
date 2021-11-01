qda_samp=function(true_status,samp,partition){
  
  n=length(true_status); p=ncol(samp)  
  error=numeric(length=partition)

  case=which(true_status==1); control=which(true_status==0)
  
  for(i in 1:partition){
    
    case_train=sort(sample(case,size=49,replace=FALSE))
    case_test=sort(setdiff(case,case_train))
    
    control_train=sort(sample(control,size=16,replace=FALSE))
    control_test=sort(setdiff(control,control_train))
    
    train=sort(c(case_train,control_train))
    test=sort(c(case_test,control_test))
    
    status_train=true_status[train]
    status_test=true_status[test]
  
    samp_test=samp[test,]
    
    s1_train=samp[case_train,]; s0_train=samp[control_train,]
    
    mu1=numeric(length=p)
    for(j in 1:49){
      mu1=mu1+as.numeric(s1_train[j,])
    }
    mu1=mu1/49
    
    mu0=numeric(length=p)
    for(l in 1:16){
      mu0=mu0+as.numeric(s0_train[l,])
    }
    mu0=mu0/16
    
    est_cov0=1/16*t(s0_train)%*%(diag(16)-1/16*rep(1,16)%*%t(rep(1,16)))%*%s0_train+0.0000001*diag(p)
    est_prec0=solve(est_cov0)
    
    est_cov1=1/49*t(s1_train)%*%(diag(49)-1/49*rep(1,49)%*%t(rep(1,49)))%*%s1_train+0.0000001*diag(p)
    est_prec1=solve(est_cov1)
    
    p1=49/65
    p0=16/65
    
    error_temp=0
    for(ind in 1:130){
    est=0
    ob=samp_test[ind,]
    rule1=as.numeric(-1/2*log(det(est_cov1))-1/2*t(ob-mu1)%*%est_prec1%*%(ob-mu1)+log(p1))
    rule0=as.numeric(-1/2*log(det(est_cov0))-1/2*t(ob-mu0)%*%est_prec0%*%(ob-mu0)+log(p0))
    
    if(rule1>rule0){
      est=1
    }
    diff=est-status_test[ind]
    error_temp=error_temp+abs(diff)
    }
  error[i]=error_temp/130
  }
  
  return(error)
}
  
 

qda_gl=function(true_status,samp,partition){
  
  n=length(true_status); p=ncol(samp)  
  error=numeric(length=partition)
  
  case=which(true_status==1); control=which(true_status==0)
  
  step.size=100
  tol=0.001
  P=matrix(1,p,p)
  diag(P)=0
  lam=0.02

  for(i in 1:partition){
    
    case_train=sort(sample(case,size=49,replace=FALSE))
    case_test=sort(setdiff(case,case_train))
    
    control_train=sort(sample(control,size=16,replace=FALSE))
    control_test=sort(setdiff(control,control_train))
    
    train=sort(c(case_train,control_train))
    test=sort(c(case_test,control_test))
    
    status_train=true_status[train]
    status_test=true_status[test]
    
    samp_test=samp[test,]
    
    s1_train=samp[case_train,]; s0_train=samp[control_train,]
    
    mu1=numeric(length=p)
    for(j in 1:49){
      mu1=mu1+as.numeric(s1_train[j,])
    }
    mu1=mu1/49
    
    mu0=numeric(length=p)
    for(l in 1:16){
      mu0=mu0+as.numeric(s0_train[l,])
    }
    mu0=mu0/16
    
    samp_cov0=1/16*t(s0_train)%*%(diag(16)-1/16*rep(1,16)%*%t(rep(1,16)))%*%s0_train+0.0001*diag(p)
    est_cov0=spcov(Sigma=samp_cov0,S=samp_cov0,lambda=lam*P,step.size=step.size,n.inner.steps=200,
                   thr.inner=0,tol.outer=tol,trace=1)$Sigma
    est_prec0=solve(est_cov0)
    
    samp_cov1=1/49*t(s1_train)%*%(diag(49)-1/49*rep(1,49)%*%t(rep(1,49)))%*%s1_train+0.0001*diag(p)
    est_cov1=spcov(Sigma=samp_cov1,S=samp_cov1,lambda=lam*P,step.size=step.size,n.inner.steps=200,
                    thr.inner=0,tol.outer=tol,trace=1)$Sigma
    est_prec1=solve(est_cov1)
    
    p1=49/65; p0=16/65
    
    error_temp=0
    for(ind in 1:130){
      est=0
      ob=samp_test[ind,]
      rule1=as.numeric(-1/2*log(det(est_cov1))-1/2*t(ob-mu1)%*%est_prec1%*%(ob-mu1)+log(p1))
      rule0=as.numeric(-1/2*log(det(est_cov0))-1/2*t(ob-mu0)%*%est_prec0%*%(ob-mu0)+log(p0))
      
      if(rule1>rule0){
        est=1
      }
      diff=est-status_test[ind]
      error_temp=error_temp+abs(diff)
    }
    error[i]=error_temp/130
  }
  return(error)
}

qda_proposed=function(true_status,samp,partition){
  
  n=length(true_status); p=ncol(samp)  
  error_mpp=numeric(length=partition)
  error_map=numeric(length=partition)
  
  case=which(true_status==1); control=which(true_status==0)
  
  for(i in 1:partition){
    
    print(paste(i,"th replication",sep=""))
    
    case_train=sort(sample(case,size=49,replace=FALSE))
    case_test=sort(setdiff(case,case_train))
    
    control_train=sort(sample(control,size=16,replace=FALSE))
    control_test=sort(setdiff(control,control_train))
    
    train=sort(c(case_train,control_train))
    test=sort(c(case_test,control_test))
    
    status_train=true_status[train]
    status_test=true_status[test]
    
    samp_test=samp[test,]
    
    s1_train=samp[case_train,]; s0_train=samp[control_train,]
    
    mu1=numeric(length=p)
    for(j in 1:49){
      mu1=mu1+as.numeric(s1_train[j,])
    }
    mu1=mu1/49
    
    mu0=numeric(length=p)
    for(l in 1:16){
      mu0=mu0+as.numeric(s0_train[l,])
    }
    mu0=mu0/16
    
    p1=49/65; p0=16/65

    samp_cov0=1/16*t(s0_train)%*%(diag(16)-1/16*rep(1,16)%*%t(rep(1,16)))%*%s0_train+0.001*diag(p)
    mcmc_samp0=mcmc(16,samp_cov0,5,16*0.06,0.1^3,12000,3000,0.6)
    mcmc_mpp0=mpp(mcmc_samp0)
    mcmc_map0=mmp(mcmc_samp0)
    
    est_cov_mpp0=bcd_covariance2(samp_cov0,16,5,16*0.06,10000,0.1^3,mcmc_mpp0)
    est_prec_mpp0=solve(est_cov_mpp0)
    
    est_cov_map0=bcd_covariance2(samp_cov0,16,5,16*0.06,10000,0.1^3,mcmc_map0)
    est_prec_map0=solve(est_cov_map0)
    
    samp_cov1=1/49*t(s1_train)%*%(diag(49)-1/49*rep(1,49)%*%t(rep(1,49)))%*%s1_train+0.001*diag(p)
    mcmc_samp1=mcmc(49,samp_cov1,5,49*0.06,0.1^3,12000,3000,0.6)
    mcmc_mpp1=mpp(mcmc_samp1)
    mcmc_map1=mmp(mcmc_samp1)
    
    est_cov_mpp1=bcd_covariance2(samp_cov1,49,1,49*0.06,10000,0.1^3,mcmc_mpp1)
    est_prec_mpp1=solve(est_cov_mpp1)
    
    est_cov_map1=bcd_covariance2(samp_cov1,49,1,49*0.06,10000,0.1^3,mcmc_map1)
    est_prec_map1=solve(est_cov_map1)
    
    error_mpp_temp=0; error_map_temp=0
    
    for(ind in 1:130){
      est_mpp=0; est_map=0
    
      ob=samp_test[ind,]
      rule1_mpp=as.numeric(-1/2*log(det(est_cov_mpp1))-1/2*t(ob-mu1)%*%est_prec_mpp1%*%(ob-mu1)+log(p1))
      rule0_mpp=as.numeric(-1/2*log(det(est_cov_mpp0))-1/2*t(ob-mu0)%*%est_prec_mpp0%*%(ob-mu0)+log(p0))
      
      if(rule1_mpp>rule0_mpp){
        est_mpp=1
      }
      
      rule1_map=as.numeric(-1/2*log(det(est_cov_map1))-1/2*t(ob-mu1)%*%est_prec_map1%*%(ob-mu1)+log(p1))
      rule0_map=as.numeric(-1/2*log(det(est_cov_map0))-1/2*t(ob-mu0)%*%est_prec_map0%*%(ob-mu0)+log(p0))
      if(rule1_map>rule0_map){
        est_map=1
      }
      diff_mpp=est_mpp-status_test[ind]
      error_mpp_temp=error_mpp_temp+abs(diff_mpp)
      
      diff_map=est_map-status_test[ind]
      error_map_temp=error_map_temp+abs(diff_map)
    }
    error_mpp_temp/130; error_map_temp/130
    print(paste(error_mpp_temp/130,error_map_temp/130,sep=" "))
    error_mpp[i]=error_mpp_temp/130
    error_map[i]=error_map_temp/130
  }
  return(list(error_mpp,error_map))
}