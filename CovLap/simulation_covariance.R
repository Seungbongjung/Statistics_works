#algorithm

bcd_concentration=function(samp,n,v0,v1,lambda,tol_it,tol,edge){
  
  p=nrow(samp)
  
  #edge penalty
  edge_penalty=numeric(length=choose(p,2))
  edge_penalty[which(edge==0)]=1/(n*v0^2)
  edge_penalty[which(edge==1)]=1/(n*v1^2)
  Lambda=matrix(rep(0,p^2),nrow=p,ncol=p)
  Lambda[upper.tri(Lambda)]=edge_penalty
  Lambda=Lambda+t(Lambda)
  
  init_omega=solve(samp)
  
  # permutation matrix for looping through columns & rows
  perms=matrix(NA,nrow=p-1,ncol=p)
  permInt=1:p
  for(i in 1:ncol(perms))
  {
    perms[,i]<-permInt[-i]
  }
  
  omega=init_omega
  int=0
  
  for(i in 1:tol_it){
    
    omega_temp=omega 
    
    for(j in 1:p){
      omega11=omega_temp[perms[,j],perms[,j]]
      s12=samp[perms[,j],j]
      s22=samp[j,j]
      Lambda12=diag(Lambda[perms[,j],j])
      
      pen=s22+lambda/n
      omega12=-solve(Lambda12+pen*solve(omega11))%*%s12
      omega22=1/pen+t(s12)%*%solve(omega11)%*%s12
      
      omega_temp[j,j]=omega22
      omega_temp[perms[,j],j]=omega12
      omega_temp[j,perms[,j]]=omega12
    }
    
    error=sqrt(sum((omega_temp-omega)^2))
    int=int+1
    if(error<tol){
      return(omega_temp)
    }else{
      omega=omega_temp
    }
  }
  return(omega)
}

# edge count
edge_ind=function(m,tol){
  edge_list=m[upper.tri(m)]
  edge_indicator=as.numeric((abs(edge_list)>=tol))
  return(edge_indicator)
}

# performance check
performance=function(esti,m,tol){
  edge_list_est=esti[upper.tri(esti)]
  edge_list_m=m[upper.tri(m)]
  
  edge_ind=as.numeric(abs(edge_list_est)>=tol)
  true_edge=as.numeric(abs(edge_list_m)>=tol)
  
  edge_ind1=which(edge_ind==1)
  edge_ind0=which(edge_ind==0)
  
  true_edge1=which(true_edge==1)
  true_edge0=which(true_edge==0)
  
  
  tp=length(intersect(edge_ind1,true_edge1))
  tn=length(intersect(edge_ind0,true_edge0))
  fp=length(intersect(edge_ind1,true_edge0))
  fn=length(intersect(edge_ind0,true_edge1))
  
  sp=tn/(tn+fp)
  se=tp/(tp+fn)
  return(c(sp,se))
}

#simulation

library(mvtnorm)
p1=30; p2=50; p3=100; n1=100; n2=200

#models used for simulation

#Model1-AR(1)
sigma=matrix(nrow=p3,ncol=p3)
for(i in 1:p3){
  for(j in 1:p3){
    sigma[i,j]=0.7^{abs(i-j)}
  }
}
model1_sigma1=sigma[1:p1,1:p1]
model1_sigma2=sigma[1:p2,1:p2]
model1_sigma3=sigma[1:p3,1:p3]

model1_sigma1_edge=edge_ind(solve(model1_sigma1),0.1^3)
model1_sigma2_edge=edge_ind(solve(model1_sigma2),0.1^3)
model1_sigma3_edge=edge_ind(solve(model1_sigma3),0.1^3)

#Model2-AR(2)
omega=diag(p3)
for(i in 1:p3){
  if(i>1){
    omega[i,i-1]=0.5
    omega[i-1,i]=0.5
  }
  if(i>2){
    omega[i,i-2]=0.25
    omega[i-2,i]=0.25
  }
}

model2_sigma1=solve(omega[1:p1,1:p1])
model2_sigma2=solve(omega[1:p2,1:p2])
model2_sigma3=solve(omega[1:p3,1:p3])

model2_sigma1_edge=edge_ind(solve(model2_sigma1),0.1^3)
model2_sigma2_edge=edge_ind(solve(model2_sigma2),0.1^3)
model2_sigma3_edge=edge_ind(solve(model2_sigma3),0.1^3)

#Model3-Star
omega_star=diag(p3)
for(i in 2:p3){
  omega_star[1,i]=0.1
  omega_star[i,1]=0.1
}

model3_sigma1=solve(omega_star[1:p1,1:p1])
model3_sigma2=solve(omega_star[1:p2,1:p2])
model3_sigma3=solve(omega_star[1:p3,1:p3])

model3_sigma1_edge=edge_ind(solve(model3_sigma1),0.1^3)
model3_sigma2_edge=edge_ind(solve(model3_sigma2),0.1^3)
model3_sigma3_edge=edge_ind(solve(model3_sigma3),0.1^3)

#Model4-Circle
omega_circle1=2*diag(p1); omega_circle2=2*diag(p2); omega_circle3=2*diag(p3)

for(i in 2:p1){
  omega_circle1[i,i-1]=1
  omega_circle1[i-1,i]=1
}
omega_circle1[1,p1]=0.9
omega_circle1[p1,1]=0.9

for(i in 2:p2){
  omega_circle2[i,i-1]=1
  omega_circle2[i-1,i]=1
}
omega_circle2[1,p2]=0.9
omega_circle2[p2,1]=0.9

for(i in 2:p3){
  omega_circle3[i,i-1]=1
  omega_circle3[i-1,i]=1
}
omega_circle3[1,p3]=0.9
omega_circle3[p3,1]=0.9

model4_sigma1=solve(omega_circle1)
model4_sigma2=solve(omega_circle2)
model4_sigma3=solve(omega_circle3)

model4_sigma1_edge=edge_ind(solve(model4_sigma1),0.1^3)
model4_sigma2_edge=edge_ind(solve(model4_sigma2),0.1^3)
model4_sigma3_edge=edge_ind(solve(model4_sigma3),0.1^3)