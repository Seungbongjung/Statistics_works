bcd_covariance1=function(samp,n,v0,v1,lambda,tol_it,tol,edge){
  
  p=nrow(samp)
  
  #edge penalty
  edge_penalty=numeric(length=choose(p,2))
  edge_penalty[which(edge==0)]=1/(n*v0^2)
  edge_penalty[which(edge==1)]=1/(n*v1^2)
  Lambda=matrix(rep(0,p^2),nrow=p,ncol=p)
  Lambda[upper.tri(Lambda)]=edge_penalty
  Lambda=Lambda+t(Lambda)
  
  rho=lambda/n
  init_cov=diag(diag(samp))+rho*diag(p)
  
  # permutation matrix for looping through columns & rows
  perms=matrix(NA,nrow=p-1,ncol=p)
  permInt=1:p
  for(i in 1:ncol(perms))
  {
    perms[,i]<-permInt[-i]
  }
  
  cov=init_cov
  k=0
  
  for(i in 1:tol_it){
    
    cov_temp=cov 
    
    for(j in 1:p){
      
      cov11=cov_temp[perms[,j],perms[,j]]; omega11=solve(cov11)
      cov12=cov_temp[perms[,j],j]
      
      s11=samp[perms[,j],perms[,j]]
      s12=samp[perms[,j],j]
      s22=samp[j,j]
    
      Lambda12=diag(Lambda[perms[,j],j])
    
      omega_s11=omega11%*%s11%*%omega11
      u=t(cov12)%*%omega_s11%*%cov12-2*t(s12)%*%omega11%*%cov12+s22
      u=as.numeric(u)
      
      gamma=(-1+sqrt(1+4*u*rho))/(2*rho)
      beta=solve(Lambda12+rho*omega11+omega_s11/gamma)%*%omega11%*%s12/gamma
      
      cov_temp[j,j]=gamma+t(beta)%*%omega11%*%beta
      cov_temp[perms[,j],j]=beta
      cov_temp[j,perms[,j]]=beta
      
    }
    
    error=sqrt(sum((cov_temp-cov)^2))
    k=k+1
    if(error<tol){
      return(cov_temp)
    }else{
      cov=cov_temp
    }
  }
  return(cov)
}

bcd_covariance2=function(samp,n,v,lambda,tol_it,tol,edge){
  
  p=nrow(samp)

  #edge penalty
  edge_penalty=rep(0,choose(p,2))
  edge_penalty[which(edge==1)]=1/(n*v^2)
  Lambda=matrix(rep(0,p^2),nrow=p,ncol=p)
  Lambda[upper.tri(Lambda)]=edge_penalty
  Lambda=Lambda+t(Lambda)
  edge_list=matrix(rep(0,p^2),nrow=p,ncol=p)
  edge_list[upper.tri(edge_list)]=edge
  edge_list=edge_list+t(edge_list)
  
  rho=lambda/n
  init_cov=diag(diag(samp))+rho*diag(p)
  
  # permutation matrix for looping through columns & rows
  perms=matrix(NA,nrow=p-1,ncol=p)
  permInt=1:p
  for(i in 1:ncol(perms)){
    perms[,i]<-permInt[-i]
  }
  
  cov=init_cov
  k=0
  
  for(i in 1:tol_it){
  
    cov_temp=cov 
    
    for(j in 1:p){
      
      cov11=cov_temp[perms[,j],perms[,j]]; omega11=solve(cov11)
      cov12=cov_temp[perms[,j],j]
      
      s11=samp[perms[,j],perms[,j]]
      s12=samp[perms[,j],j]
      s22=samp[j,j]
    
      edge12=as.vector(edge_list[perms[,j],j])
      edge12_list=which(edge12==1)
      Lambda12=diag(Lambda[perms[,j],j])
      
      omega_s11=omega11%*%s11%*%omega11
      u=t(cov12)%*%omega_s11%*%cov12-2*t(s12)%*%omega11%*%cov12+s22
      u=as.numeric(u)
      
      gamma=(-1+sqrt(1+4*u*rho))/(2*rho)
      
      if(length(edge12_list)!=0){
        quad_mat=(Lambda12+rho*omega11+omega_s11/gamma)[edge12_list,edge12_list]
        obj_vec=(omega11%*%s12)[edge12_list]
        beta=solve(quad_mat)%*%obj_vec/gamma
        cov12[edge12_list]=beta
      }
    
      cov_temp[j,j]=gamma+t(cov12)%*%omega11%*%cov12
      cov_temp[perms[,j],j]=cov12
      cov_temp[j,perms[,j]]=cov12
      
    }
    error=sqrt(sum((cov_temp-cov)^2))
    k=k+1
    if(error<tol){
      return(cov)
    }else{
      cov=cov_temp
    }
  }
  return(cov)
}
