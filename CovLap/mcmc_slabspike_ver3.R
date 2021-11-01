source("bcd_algorithm.R")
source("performance.R")

objec=function(m,n,samp,edge,v,lambda){
  m_upper=m[upper.tri(m)]
  ind=which(edge==1)
  if(length(ind)>0){
    m_upper=m_upper[which(edge==1)] 
  }else{
    m_upper=0
  }
  omega=solve(m)
  value=log(det(m))+sum(diag(samp%*%omega))+1/(n*v^2)*sum(m_upper^2)+lambda/n*sum(diag(m))
  return(value)
}

hessian=function(cov,s,edge,n,v){
  ind=which(edge==1)
  edge_no=sum(edge)
  p=nrow(cov)
  omega=solve(cov)
  u=omega%*%s%*%omega
  
  E=array(NA,dim=c(1,2,0))
  for(k in 1:(p-1)){
    for(l in (k+1):p){
      E=array(c(E,c(l,k)),dim=dim(E)+c(0,0,1))
    }
  }
  E=matrix(E,ncol=2,byrow=TRUE)
  edge_list=E[ind,]
  
  h=matrix(nrow=p+edge_no,ncol=p+edge_no)
  
  if(edge_no==0){
    for(i in 1:p){
      for(l in i:p){
        h[i,l]=2*u[i,l]*omega[l,i]-(omega[i,l])^2
      }
    }
    h=h+t(h)-diag(diag(h))
    return(h)
  }else{
    h11=matrix(rep(0,p^2),nrow=p)
    for(i in 1:p){
      for(l in i:p){
        h11[i,l]=2*u[i,l]*omega[l,i]-(omega[i,l])^2
      }
    }
    h11=h11+t(h11)-diag(diag(h11))
    h[1:p,1:p]=h11
    
    h12=matrix(rep(0,p*edge_no),nrow=p)
    for(i in 1:p){
      for(l in 1:edge_no){
        a=edge_list[l,1]; b=edge_list[l,2]
        h12[i,l]=2*(-omega[a,i]*omega[b,i]+u[a,i]*omega[b,i]+u[b,i]*omega[a,i])
      }
    }
    h[1:p,(p+1):(p+edge_no)]=h12
    h[(p+1):(p+edge_no),1:p]=t(h12)
    
    h22=matrix(rep(0,edge_no^2),nrow=edge_no)
    for(i in 1:edge_no){
      for(l in i:edge_no){
        a=edge_list[i,1]; b=edge_list[i,2]
        c=edge_list[l,1]; d=edge_list[l,2]
        h22[i,l]=2*(-omega[a,c]*omega[b,d]-omega[a,d]*omega[b,c]
                    +u[b,c]*omega[a,d]+u[b,d]*omega[a,c]+u[c,a]*omega[b,d]+u[d,a]*omega[b,c])
      }
    }
    h22=h22+t(h22)-diag(diag(h22))+diag(rep(2/(n*v^2),edge_no))
    h[(p+1):(p+edge_no),(p+1):(p+edge_no)]=h22
    return(h)
  }
}

post_ratio=function(m1,m2,n,samp,edge1,edge2,v,lambda,q){
  post1=objec(m1,n,samp,edge1,v,lambda); post2=objec(m2,n,samp,edge2,v,lambda)
  hessian1=hessian(m1,s,edge1,n,v); hessian2=hessian(m2,s,edge2,n,v)
  
  post_diff=post1-post2
  edge_diff=sum(edge1)-sum(edge2)
  
  ratio=det(hessian1)/det(hessian2)
  
  result=(q/((1-q)*v*sqrt(2*pi))*sqrt(4*pi/n))^edge_diff*exp(-n/2*post_diff)*1/sqrt(ratio)
  return(result)
}

unif_samp=function(edge,p){
  l=choose(p,2)
  samp=sample(1:l,1)
  result=edge
  if(edge[samp]==1){
    result[samp]=0
  }
  if(edge[samp]==0){
    result[samp]=1
  }
  return(result)
}

mcmc=function(n,s,v,lambda,tol,size,burnin,q){
  
  p=nrow(s)
  edge_temp=edge_ind(s,0.1^3)
  edge_temp_ind=which(edge_temp==1)
  edge_temp_ind=sample(edge_temp_ind,100,replace=FALSE)
  edge_temp[-edge_temp_ind]=0
  cov_temp=bcd_covariance2(s,n,v,lambda,10000,0.1^3,edge_temp)
  
  total_samp=matrix(rep(0,size*choose(p,2)),nrow=size)
  
  for(i in 1:(size+burnin)){
    if(i<=burnin){
      edge_cand=unif_samp(edge_temp,p)
      cov_cand=bcd_covariance2(s,n,v,lambda,10000,tol,edge_cand)
      alpha=min(c(1,post_ratio(cov_cand,cov_temp,n,s,edge_cand,edge_temp,v,lambda,q)))
      u=runif(1,0,1)
    
      if(u<=alpha){
        edge_temp=edge_cand
        cov_temp=cov_cand
      }
      if(i%%100==0){
        print(paste(i,"th burnin",sep=""))
      }
    }else{
      edge_cand=unif_samp(edge_temp,p)
      cov_cand=bcd_covariance2(s,n,v,lambda,10000,tol,edge_cand)
      alpha=min(c(1,post_ratio(cov_cand,cov_temp,n,s,edge_cand,edge_temp,v,lambda,q)))
      u=runif(1,0,1)
      if(u<=alpha){
        edge_temp=edge_cand
        cov_temp=cov_cand
      }
      total_samp[i-burnin,]=edge_temp
      if(i%%100==0){
        print(paste(i-burnin,"th sample",sep=""))
      }
    }
  }
  return(total_samp)
}

