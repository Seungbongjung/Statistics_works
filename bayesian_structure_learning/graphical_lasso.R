#simulation

library(tidyverse)
library(mvtnorm)
library(glasso)

#dimension and sample size
p1=30; p2=50; p3=100
n1=100; n2=200

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

#Model3-Star
omega_star=diag(p3)
for(i in 2:p3){
  omega_star[1,i]=0.1
  omega_star[i,1]=0.1
}

model3_sigma1=solve(omega_star[1:p1,1:p1])
model3_sigma2=solve(omega_star[1:p2,1:p2])
model3_sigma3=solve(omega_star[1:p3,1:p3])

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

#Edge count

edge_count=function(m){
  temp=m[upper.tri(m)]
  temp=abs(temp)>0.1^3
  return(temp)
}

#Simulation of graphical lasso
graphical_lasso=function(n,p,sigma,rho,iter,seed){
  edge_temp=edge_count(solve(sigma))
  result=list(tp=numeric(iter),tn=numeric(iter),fp=numeric(iter),fn=numeric(iter))
  set.seed(seed)
  for(i in 1:iter){
    samp=rmvnorm(n,mean=rep(0,p),sigma=sigma)
    s=t(samp)%*%samp/n
    gl=glasso(s=s,rho=rho,approx=FALSE)
    estimated_var=gl$wi
    estimated_var_temp=edge_count(estimated_var)
    result[[1]][i]=length(intersect(which(estimated_var_temp==1),which(edge_temp==1)))
    result[[2]][i]=length(intersect(which(estimated_var_temp==0),which(edge_temp==0)))
    result[[3]][i]=length(intersect(which(estimated_var_temp==1),which(edge_temp==0)))
    result[[4]][i]=length(intersect(which(estimated_var_temp==0),which(edge_temp==1)))
  }
  result[[1]][i]=length(intersect(which(estimated_var_temp==1),which(edge_temp==1)))
  result[[2]][i]=length(intersect(which(estimated_var_temp==0),which(edge_temp==0)))
  result[[3]][i]=length(intersect(which(estimated_var_temp==1),which(edge_temp==0)))
  result[[4]][i]=length(intersect(which(estimated_var_temp==0),which(edge_temp==1)))
  return(result)
}

#Accuracy 
accuracy=function(data){
  tp=data[[1]]; tn=data[[2]]; fp=data[[3]]; fn=data[[4]]
  sp=tn/(tn+fp)
  se=tp/(tp+fn)
  mcc=(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  fpr=1-sp
  return(list(sp,se,mcc,fpr))
}

#Model1-simulation 
model1_sigma1_samp1=graphical_lasso(n1,p1,sigma=model1_sigma1,rho=0.5,iter=100,seed=1)
model1_sigma1_samp1_acc=accuracy(model1_sigma1_samp1)

model1_sigma1_samp2=graphical_lasso(n2,p1,sigma=model1_sigma1,rho=0.5,iter=100,seed=2)
model1_sigma1_samp2_acc=accuracy(model1_sigma1_samp2)

model1_sigma2_samp1=graphical_lasso(n1,p2,sigma=model1_sigma2,rho=0.5,iter=100,seed=3)
model1_sigma2_samp1_acc=accuracy(model1_sigma2_samp1)

model1_sigma2_samp2=graphical_lasso(n2,p2,sigma=model1_sigma2,rho=0.5,iter=100,seed=4)
model1_sigma2_samp2_acc=accuracy(model1_sigma2_samp2)

model1_sigma3_samp1=graphical_lasso(n1,p3,sigma=model1_sigma3,rho=0.5,iter=100,seed=5)
model1_sigma3_samp1_acc=accuracy(model1_sigma3_samp1)

model1_sigma3_samp2=graphical_lasso(n2,p3,sigma=model1_sigma3,rho=0.5,iter=100,seed=6)
model1_sigma3_samp2_acc=accuracy(model1_sigma3_samp2)

#Model2-simulation 
model2_sigma1_samp1=graphical_lasso(n1,p1,sigma=model2_sigma1,rho=0.5,iter=100,seed=7)
model2_sigma1_samp1_acc=accuracy(model2_sigma1_samp1)

model2_sigma1_samp2=graphical_lasso(n2,p1,sigma=model2_sigma1,rho=0.5,iter=100,seed=8)
model2_sigma1_samp2_acc=accuracy(model2_sigma1_samp2)

model2_sigma2_samp1=graphical_lasso(n1,p2,sigma=model2_sigma2,rho=0.5,iter=100,seed=9)
model2_sigma2_samp1_acc=accuracy(model2_sigma2_samp1)

model2_sigma2_samp2=graphical_lasso(n2,p2,sigma=model2_sigma2,rho=0.5,iter=100,seed=10)
model2_sigma2_samp2_acc=accuracy(model2_sigma2_samp2)

model2_sigma3_samp1=graphical_lasso(n1,p3,sigma=model2_sigma3,rho=0.5,iter=100,seed=11)
model2_sigma3_samp1_acc=accuracy(model2_sigma3_samp1)

model2_sigma3_samp2=graphical_lasso(n2,p3,sigma=model2_sigma3,rho=0.5,iter=100,seed=12)
model2_sigma3_samp2_acc=accuracy(model2_sigma3_samp2)

#Model3-Simulation
model3_sigma1_samp1=graphical_lasso(n1,p1,sigma=model3_sigma1,rho=0.2,iter=100,seed=13)
model3_sigma1_samp1_acc=accuracy(model3_sigma1_samp1)

model3_sigma1_samp2=graphical_lasso(n2,p1,sigma=model3_sigma1,rho=0.2,iter=100,seed=14)
model3_sigma1_samp2_acc=accuracy(model3_sigma1_samp2)

model3_sigma2_samp1=graphical_lasso(n1,p2,sigma=model3_sigma2,rho=0.2,iter=100,seed=15)
model3_sigma2_samp1_acc=accuracy(model3_sigma2_samp1)

model3_sigma2_samp2=graphical_lasso(n2,p2,sigma=model3_sigma2,rho=0.2,iter=100,seed=16)
model3_sigma2_samp2_acc=accuracy(model3_sigma2_samp2)

model3_sigma3_samp1=graphical_lasso(n1,p3,sigma=model3_sigma3,rho=0.2,iter=100,seed=17)
model3_sigma3_samp1_acc=accuracy(model3_sigma3_samp1)

model3_sigma3_samp2=graphical_lasso(n2,p3,sigma=model3_sigma3,rho=0.2,iter=100,seed=18)
model3_sigma3_samp2_acc=accuracy(model3_sigma3_samp2)

#Model4-simulation

model4_sigma1_samp1=graphical_lasso(n1,p1,sigma=model4_sigma1,rho=0.5,iter=100,seed=19)
model4_sigma1_samp1_acc=accuracy(model4_sigma1_samp1)

model4_sigma1_samp2=graphical_lasso(n2,p1,sigma=model4_sigma1,rho=0.5,iter=100,seed=20)
model4_sigma1_samp2_acc=accuracy(model4_sigma1_samp2)

model4_sigma2_samp1=graphical_lasso(n1,p2,sigma=model4_sigma2,rho=0.5,iter=100,seed=21)
model4_sigma2_samp1_acc=accuracy(model4_sigma2_samp1)

model4_sigma2_samp2=graphical_lasso(n2,p2,sigma=model4_sigma2,rho=0.5,iter=100,seed=22)
model4_sigma2_samp2_acc=accuracy(model4_sigma2_samp2)

model4_sigma3_samp1=graphical_lasso(n1,p3,sigma=model4_sigma3,rho=0.5,iter=100,seed=23)
model4_sigma3_samp1_acc=accuracy(model4_sigma3_samp1)

model4_sigma3_samp2=graphical_lasso(n2,p3,sigma=model4_sigma3,rho=0.5,iter=100,seed=24)
model4_sigma3_samp2_acc=accuracy(model4_sigma3_samp2)

#Implementation of graphical lasso in bayesian context
#setting of indicator

c_gamma=function(n,p,gamma1,gamma2,lambda,q,r){
  gamma1=sum(gamma1)
  gamma2=sum(gamma2)
  gamma=gamma1-gamma2
  rho=lambda*n
  result=q^gamma*(1-q)^(-gamma)*(rho/2)^(gamma)*as.numeric(gamma1<=r)
  return(result)
}

h_gamma=function(s,prec,lambda,n){
  rho=lambda*n
  temp=prec[upper.tri(prec)]
  temp=abs(temp[which(abs(temp)>0.1^3)])
  result=-log(det(prec))+sum(diag(s%*%prec))+rho/n*sum(diag(prec))+2*rho/n*sum(temp)
  return(result)
}

H=function(mat,p){
  mat_inv=solve(mat)
  ind=matrix(rep(1:p,2),ncol=2)
  temp=mat
  temp[lower.tri(temp,diag=TRUE)]=0
  mat_ind=which(abs(temp)>0.1^3,arr.ind=TRUE)
  gamma_ind=rbind(ind,mat_ind)
  n=nrow(gamma_ind)
  result=matrix(rep(0,n^2),ncol=n)
  for(a in 1:n){
    for(b in a:n){
      i=gamma_ind[a,][1]; j=gamma_ind[a,][2]; l=gamma_ind[b,][1]; m=gamma_ind[b,][2]
      if((i!=j) & (l!=m)){
        result[a,b]=2*(mat_inv[i,l]*mat_inv[j,m]+mat_inv[i,m]*mat_inv[j,l])
      }
      if((i==j) & (l!=m)){
        result[a,b]=2*mat_inv[i,l]*mat_inv[i,m]
      }
      if((i==j) & (l==m)){
        result[a,b]=(mat_inv[i,l])^2
      }
    }
  }
  result[lower.tri(result)]=t(result)[lower.tri(result)]
  return(result)
}

posterior_ratio=function(n,p,s,prec1,prec2,q,lambda,gamma1,gamma2,r){
 C=c_gamma(n,p,gamma1,gamma2,lambda,q,r)
 h=h_gamma(s,prec1,lambda,n)-h_gamma(s,prec2,lambda,n)
 ind=sum(gamma1)-sum(gamma2)
 e1=eigen(H(prec2,p))$values; e2=eigen(H(prec1,p))$values
 e1=prod(sqrt(abs(e1))); e2=prod(sqrt(abs(e2)))
 H=e1/e2
 result=C*(4*pi/n)^(ind/2)*H*exp(-n*h/2)
 return(result)
}

q=function(gamma,p){
  m=choose(p,2)
  samp=sample(1:m,1)
  result=gamma
  if(gamma[samp]==1){
    result[samp]=0
  }
  if(gamma[samp]==0){
    result[samp]=1
  }
  return(result)
}

bayesian_GL=function(n,p,sigma,q,lambda,seed,iter1,iter2){
  set.seed(seed)
  m1=array(NA,dim=c(1,choose(p,2),iter1))
  E_p=array(NA,dim=c(1,2,0))
  for(k in 2:p){
    for(l in 1:(k-1)){
      E_p=array(c(E_p,c(l,k)),dim=dim(E_p)+c(0,0,1))
    }
  }
  E_p=matrix(E_p,ncol=2,byrow=TRUE)
  for(i in 1:iter1){
   samp=rmvnorm(n,mean=rep(0,p),sigma=sigma)
   s=1/n*t(samp)%*%samp
   prec=glasso(s,rho=lambda)$wi
   gamma=edge_count(prec)
   gamma_int=gamma
   r=sum(gamma)
   m2=array(NA,dim=c(1,choose(p,2),iter2))
   for(j in 1:iter2){
     gamma_temp=q(gamma,p)
     prec_temp=glasso(s,rho=lambda,zero=E_p[as.logical(1-gamma_temp),])$wi
     prec_gamma=edge_count(prec_temp)
     u=runif(1,0,1)
     alpha=posterior_ratio(n,p,s,prec_temp,prec,q,lambda,gamma_temp,gamma,r)
     alpha=min(c(1,alpha))
     if(u<=alpha){
       gamma=gamma_temp
       prec=prec_temp
     }
     m2[,,j]=gamma
   }
   z=iter2-100
   m2_temp=m2[,,101:iter2]
   m_temp=rowSums(m2_temp,dims=1)
   m_temp=as.vector(m_temp)/z
   m_temp=as.numeric(m_temp>=0.5)
   m1[,,i]=m_temp
  }
  return(m1)
}

model1_sigma1_samp1_bayes=bayesian_GL(n1,p1,model1_sigma1,0.4,0.5,1,50,500)
model1_sigma1_samp2_bayes=bayesian_GL(n2,p1,model1_sigma1,0.4,0.5,2,50,500)
model1_sigma2_samp1_bayes=bayesian_GL(n1,p2,model1_sigma2,0.4,0.5,3,50,500)
model1_sigma2_samp2_bayes=bayesian_GL(n2,p2,model1_sigma2,0.4,0.5,4,50,500)
model1_sigma3_samp1_bayes=bayesian_GL(n1,p3,model1_sigma3,0.4,0.5,5,50,500)
model1_sigma3_samp2_bayes=bayesian_GL(n2,p3,model1_sigma3,0.4,0.5,6,50,500)
model2_sigma1_samp1_bayes=bayesian_GL(n1,p1,model2_sigma1,0.4,0.5,7,50,500)
model2_sigma1_samp2_bayes=bayesian_GL(n2,p1,model2_sigma1,0.4,0.5,8,50,500)
model2_sigma2_samp1_bayes=bayesian_GL(n1,p2,model2_sigma2,0.4,0.5,9,50,500)
model2_sigma2_samp2_bayes=bayesian_GL(n2,p2,model2_sigma2,0.4,0.5,10,50,500)
model2_sigma3_samp1_bayes=bayesian_GL(n1,p3,model2_sigma3,0.4,0.5,11,50,500)
model2_sigma3_samp2_bayes=bayesian_GL(n2,p3,model2_sigma3,0.4,0.5,12,50,500)

model1_sigma1_samp1_bayes=matrix(model1_sigma1_samp1_bayes,ncol=50)
model1_sigma1_samp1_bayes=t(model1_sigma1_samp1_bayes)
write.csv(model1_sigma1_samp1_bayes,"model1_sigma1_samp1_bayes.csv")

model1_sigma1_samp2_bayes=matrix(model1_sigma1_samp2_bayes,ncol=50)
model1_sigma1_samp2_bayes=t(model1_sigma1_samp2_bayes)
write.csv(model1_sigma1_samp2_bayes,"model1_sigma1_samp2_bayes.csv")

model1_sigma2_samp1_bayes=matrix(model1_sigma2_samp1_bayes,ncol=50)
model1_sigma2_samp1_bayes=t(model1_sigma2_samp1_bayes)
write.csv(model1_sigma2_samp1_bayes,"model1_sigma2_samp1_bayes.csv")

model1_sigma2_samp2_bayes=matrix(model1_sigma2_samp2_bayes,ncol=50)
model1_sigma2_samp2_bayes=t(model1_sigma2_samp2_bayes)
write.csv(model1_sigma2_samp2_bayes,"model1_sigma2_samp2_bayes.csv")

model2_sigma1_samp1_bayes=matrix(model2_sigma1_samp1_bayes,ncol=50)
model2_sigma1_samp1_bayes=t(model2_sigma1_samp1_bayes)
write.csv(model2_sigma1_samp1_bayes,"model2_sigma1_samp1_bayes.csv")

model2_sigma1_samp2_bayes=matrix(model2_sigma1_samp2_bayes,ncol=50)
model2_sigma1_samp2_bayes=t(model2_sigma1_samp2_bayes)
write.csv(model2_sigma1_samp2_bayes,"model2_sigma1_samp2_bayes.csv")

model2_sigma2_samp1_bayes=matrix(model2_sigma2_samp1_bayes,ncol=50)
model2_sigma2_samp1_bayes=t(model2_sigma2_samp1_bayes)
write.csv(model2_sigma2_samp1_bayes,"model2_sigma2_samp1_bayes.csv")

model2_sigma2_samp2_bayes=matrix(model2_sigma2_samp2_bayes,ncol=50)
model2_sigma2_samp2_bayes=t(model2_sigma2_samp2_bayes)
write.csv(model2_sigma2_samp2_bayes,"model2_sigma2_samp2_bayes.csv")

model4_sigma2_samp2_bayes=bayesian_GL(n2,p2,model4_sigma2,0.4,0.5,22,50,500)
install.packages("sparsebn")
library(sparsebn)