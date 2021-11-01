lda_samp=function(true_status,samp,partition){
  
  n=length(true_status); p=ncol(samp)  
  error=numeric(length=partition)

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
    
    est_cov=1/190*t(samp_train)%*%(diag(190)-1/190*rep(1,190)%*%t(rep(1,190)))%*%samp_train+0.1^3*diag(p)
    est_prec=solve(est_cov)
    
    ob=samp[ind,]
    rule1=as.numeric(t(ob)%*%est_prec%*%mu1-1/2*t(mu1)%*%est_prec%*%mu1+log(p1))
    rule0=as.numeric(t(ob)%*%est_prec%*%mu0-1/2*t(mu0)%*%est_prec%*%mu0+log(p0))
    
    p1=71/190
    p0=119/190
    
    error_temp=0
    for(ind in 1:378){
    est=0
    ob=samp_test[ind,]
    rule1=as.numeric(t(ob)%*%est_prec%*%mu1-1/2*t(mu1)%*%est_prec%*%mu1+log(p1))
    rule0=as.numeric(t(ob)%*%est_prec%*%mu0-1/2*t(mu0)%*%est_prec%*%mu0+log(p0))
    
    if(rule1>rule0){
      est=1
    }
    diff=est-status_test[ind]
    error_temp=error_temp+abs(diff)
    }
  error[i]=error_temp/378
  }
  return(error)
}
  
 

lda_gl=function(true_status,samp,partition){
  
  n=length(true_status); p=ncol(samp)  
  error=numeric(length=partition)
  
  case=which(true_status==1); control=which(true_status==0)
  
  step.size=100
  tol=0.002
  P=matrix(1,p,p)
  diag(P)=0
  lam=0.06

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
    
    samp_cov=1/190*t(samp_train)%*%(diag(190)-1/190*rep(1,190)%*%t(rep(1,190)))%*%samp_train+0.001*diag(p)
    est_cov=spcov(Sigma=samp_cov,S=samp_cov,lambda=lam*P,step.size=step.size,n.inner.steps=200,
                   thr.inner=0,tol.outer=tol,trace=1)$Sigma
    est_prec=solve(est_cov)
    
    p1=71/190; p0=119/190
    
    error_temp=0
    for(ind in 1:378){
      est=0
      ob=samp_test[ind,]
      rule1=as.numeric(t(ob)%*%est_prec%*%mu1-1/2*t(mu1)%*%est_prec%*%mu1+log(p1))
      rule0=as.numeric(t(ob)%*%est_prec%*%mu0-1/2*t(mu0)%*%est_prec%*%mu0+log(p0))
      
      if(rule1>rule0){
        est=1
      }
      diff=est-status_test[ind]
      error_temp=error_temp+abs(diff)
    }
    error[i]=error_temp/378
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
    
    p1=71/190; p0=119/190
    
    samp_cov=1/190*t(samp_train)%*%(diag(190)-1/190*rep(1,190)%*%t(rep(1,190)))%*%samp_train+0.0001*diag(p)
    mcmc_samp=mcmc(190,samp_cov,1,190*0.05,0.1,12000,3000,0.65)
    mcmc_mpp=mpp(mcmc_samp)
    mcmc_map=mmp(mcmc_samp)

    est_cov_mpp=bcd_covariance2(samp_cov,190,1,190*0.05,10000,0.1^3,mcmc_mpp)
    est_prec_mpp=solve(est_cov_mpp)
    
    est_cov_map=bcd_covariance2(samp_cov,190,1,190*0.05,10000,0.1^3,mcmc_map)
    est_prec_map=solve(est_cov_map)
    
    error_mpp_temp=0; error_map_temp=0
    
    for(ind in 1:378){
      est_mpp=0; est_map=0
    
      ob=samp_test[ind,]
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
      diff_mpp=est_mpp-status_test[ind]
      error_mpp_temp=error_mpp_temp+abs(diff_mpp)
      
      diff_map=est_map-status_test[ind]
      error_map_temp=error_map_temp+abs(diff_map)
    }
    error_mpp_temp/378; error_map_temp/378
    print(paste(error_mpp_temp/378,error_map_temp/378,sep=" "))
    error_mpp[i]=error_mpp_temp/378
    error_map[i]=error_map_temp/378
  }
  return(list(error_mpp,error_map))
}