library(mvtnorm)
library(glasso)

p1=30; p2=50; p3=100; n1=100; n2=200

#models used for simulation

#Model1-AR(1)
sigma=matrix(nrow=p3,ncol=p3)
for(i in 1:p3){
  for(j in 1:p3){
    sigma[i,j]=0.7^{abs(i-j)}
  }
}
qnorm(0.99)
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

#calculation of posterior

c_gamma=function(n,p,gamma1,gamma2,lambda,q){
  gamma1=sum(gamma1)
  gamma2=sum(gamma2)
  gamma=gamma1-gamma2
  rho=lambda*n
  result=q^gamma*(1-q)^(-gamma)*(rho/2)^(gamma)
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

posterior_ratio=function(n,p,s,prec1,prec2,q,lambda,gamma1,gamma2){
  C=c_gamma(n,p,gamma1,gamma2,lambda,q)
  h=h_gamma(s,prec1,lambda,n)-h_gamma(s,prec2,lambda,n)
  ind=sum(gamma1)-sum(gamma2)
  e1=eigen(H(prec2,p))$values; e2=eigen(H(prec1,p))$values
  e1=sqrt(abs(e1)); e2=sqrt(abs(e2))
  l1=length(e1); l2=length(e2)
  
  if(l1>l2){
    H=prod(e1[1:l2]/e2)
    H=H*prod(e1[(l2+1):l1])
  }
  if(l1==l2){
    H=prod(e1/e2)
  }
  if(l1<l2){
    H=prod(e1/e2[1:l1])
    H=H/prod(e2[(l1+1):l2])
  }
  result=C*(4*pi/n)^(ind/2)*H*exp(-n*h/2)
  return(result)
}

edge_count=function(m){
  temp=m[upper.tri(m)]
  temp=abs(temp)>0.1^3
  return(temp)
}

bayesian_GL=function(n,p,sigma,q,lambda,seed,iter1,iter2){
  
  set.seed(seed)
  E_p=array(NA,dim=c(1,2,0))
  for(k in 2:p){
    for(l in 1:(k-1)){
      E_p=array(c(E_p,c(l,k)),dim=dim(E_p)+c(0,0,1))
    }
  }
  
  #Indicator of edge
  E_p=matrix(E_p,ncol=2,byrow=TRUE)
  
  model=list()
  edge=matrix(nrow=choose(p,2),ncol=iter1)
  
  for(i in 1:iter1){
    #iter1 refers to number of replications  
    
    samp=rmvnorm(n,mean=rep(0,p),sigma=sigma)
    s=1/n*t(samp)%*%samp
    prec=glasso(s,rho=lambda)$wi
    gamma=edge_count(prec)
    r=sum(gamma)-2
    edge_list=which(gamma==1)
    edge_change=sample(edge_list,3)
    gamma_final=gamma
    gamma_final[edge_change]=0
    edge_temp=E_p[as.logical(gamma_final),]
    prec_temp=matrix(rep(0,p^2),nrow=p)
    
    l=nrow(edge_temp)
    for(u in 1:l){
      prec_temp[edge_temp[u,1],edge_temp[u,2]]=
        prec[edge_temp[u,1],edge_temp[u,2]]
    }
    
    diag(prec_temp)=diag(prec)
    prec_temp[lower.tri(prec_temp)]=
      t(prec_temp)[lower.tri(t(prec_temp))]
    prec_final=prec_temp

    for(j in 1:iter2){
      
      gamma_temp=gamma
      edge_alter=sample(1:r,1)

      edge_sub=sample(edge_list,edge_alter)
      gamma_temp[edge_sub]=0
      edge_temp=E_p[as.logical(gamma_temp),]
      prec_temp=matrix(rep(0,p^2),nrow=p)
      l=nrow(edge_temp)
      
      for(u in 1:l){
        prec_temp[edge_temp[u,1],edge_temp[u,2]]=
          prec[edge_temp[u,1],edge_temp[u,2]]
      }
      diag(prec_temp)=diag(prec)
      prec_temp[lower.tri(prec_temp)]=
        t(prec_temp)[lower.tri(t(prec_temp))]
      
      post_ratio=posterior_ratio(n,p,s,prec_temp,prec_final,q,lambda,
                                 gamma_temp,gamma_final)
      if(post_ratio>=1){
        prec_final=prec_temp
        gamma_final=gamma_temp
      }
    }
    model[[i]]=prec_final
    edge[,i]=gamma_final
    if(i%%10 == 0){
      print(paste0("Replication : ",i))
    }
  }
  return(list(model,edge))
}

#Model1-AR(1) model

model1_sigma1_bayes1=bayesian_GL(n1,p1,model1_sigma1,0.45,0.5,1,300,100)
model1_sigma1_bayes2=bayesian_GL(n2,p1,model1_sigma1,0.45,0.5,2,300,100)
model1_sigma2_bayes1=bayesian_GL(n1,p2,model1_sigma2,0.45,0.5,3,300,100)
model1_sigma2_bayes2=bayesian_GL(n2,p2,model1_sigma2,0.45,0.5,4,300,100)
model1_sigma3_bayes1=bayesian_GL(n1,p3,model1_sigma3,0.45,0.5,5,300,100)
model1_sigma3_bayes2=bayesian_GL(n2,p3,model1_sigma3,0.45,0.5,6,300,100)

write.csv(model1_sigma1_bayes1[[2]],"model1_sigma1_bayes1.csv")
write.csv(model1_sigma1_bayes2[[2]],"model1_sigma1_bayes2.csv")
write.csv(model1_sigma2_bayes1[[2]],"model1_sigma2_bayes1.csv")
write.csv(model1_sigma2_bayes2[[2]],"model1_sigma2_bayes2.csv")
write.csv(model1_sigma3_bayes1[[2]],"model1_sigma3_bayes1.csv")
write.csv(model1_sigma3_bayes2[[2]],"model1_sigma3_bayes2.csv")

#Model2-AR(2) model

model2_sigma1_bayes1=bayesian_GL(n1,p1,model2_sigma1,0.45,0.5,7,300,100)
model2_sigma1_bayes2=bayesian_GL(n2,p1,model2_sigma1,0.45,0.5,8,300,100)
model2_sigma2_bayes1=bayesian_GL(n1,p2,model2_sigma2,0.45,0.5,9,300,100)
model2_sigma2_bayes2=bayesian_GL(n2,p2,model2_sigma2,0.45,0.5,10,300,100)
model2_sigma3_bayes1=bayesian_GL(n1,p3,model2_sigma3,0.45,0.5,11,300,100)
model2_sigma3_bayes2=bayesian_GL(n2,p3,model2_sigma3,0.45,0.5,12,300,100)

write.csv(model2_sigma1_bayes1[[2]],"model2_sigma1_bayes1.csv")
write.csv(model2_sigma1_bayes2[[2]],"model2_sigma1_bayes2.csv")
write.csv(model2_sigma2_bayes1[[2]],"model2_sigma2_bayes1.csv")
write.csv(model2_sigma2_bayes2[[2]],"model2_sigma2_bayes2.csv")
write.csv(model2_sigma3_bayes1[[2]],"model2_sigma3_bayes1.csv")
write.csv(model2_sigma3_bayes2[[2]],"model2_sigma3_bayes2.csv")

#Model3-Star model

model3_sigma1_bayes1=bayesian_GL(n1,p1,model3_sigma1,0.45,0.2,13,300,100)
model3_sigma1_bayes2=bayesian_GL(n2,p1,model3_sigma1,0.45,0.2,14,300,100)
model3_sigma2_bayes1=bayesian_GL(n1,p2,model3_sigma2,0.45,0.2,15,300,100)
model3_sigma2_bayes2=bayesian_GL(n2,p2,model3_sigma2,0.45,0.2,16,300,100)
model3_sigma3_bayes1=bayesian_GL(n1,p3,model3_sigma3,0.45,0.2,17,300,100)
model3_sigma3_bayes2=bayesian_GL(n2,p3,model3_sigma3,0.45,0.2,18,300,100)

write.csv(model3_sigma1_bayes1[[2]],"model3_sigma1_bayes1.csv")
write.csv(model3_sigma1_bayes2[[2]],"model3_sigma1_bayes2.csv")
write.csv(model3_sigma2_bayes1[[2]],"model3_sigma2_bayes1.csv")
write.csv(model3_sigma2_bayes2[[2]],"model3_sigma2_bayes2.csv")
write.csv(model3_sigma3_bayes1[[2]],"model3_sigma3_bayes1.csv")
write.csv(model3_sigma3_bayes2[[2]],"model3_sigma3_bayes2.csv")

#Model4-Circle model 

model4_sigma1_bayes1=bayesian_GL(n1,p1,model4_sigma1,0.45,0.5,19,300,100)
model4_sigma1_bayes2=bayesian_GL(n2,p1,model4_sigma1,0.45,0.5,20,300,100)
model4_sigma2_bayes1=bayesian_GL(n1,p2,model4_sigma2,0.45,0.5,21,300,100)
model4_sigma2_bayes2=bayesian_GL(n2,p2,model4_sigma2,0.45,0.5,22,300,100)
model4_sigma3_bayes1=bayesian_GL(n1,p3,model4_sigma3,0.45,0.5,23,1,1)
model4_sigma3_bayes2=bayesian_GL(n2,p3,model4_sigma3,0.45,0.5,24,300,100)

write.csv(model4_sigma1_bayes1[[2]],"model4_sigma1_bayes1.csv")
write.csv(model4_sigma1_bayes2[[2]],"model4_sigma1_bayes2.csv")
write.csv(model4_sigma2_bayes1[[2]],"model4_sigma2_bayes1.csv")
write.csv(model4_sigma2_bayes2[[2]],"model4_sigma2_bayes2.csv")
write.csv(model4_sigma3_bayes1[[2]],"model4_sigma3_bayes1.csv")
write.csv(model4_sigma3_bayes2[[2]],"model4_sigma3_bayes2.csv")

### Performance

model1_sigma1_bayes1=as.matrix(read.csv("model1_sigma1_bayes1.csv")[,-1])
model1_sigma1_bayes2=as.matrix(read.csv("model1_sigma1_bayes2.csv")[,-1])
model1_sigma2_bayes1=as.matrix(read.csv("model1_sigma2_bayes1.csv")[,-1])
model1_sigma2_bayes2=as.matrix(read.csv("model1_sigma2_bayes2.csv")[,-1])
model1_sigma3_bayes1=as.matrix(read.csv("model1_sigma3_bayes1.csv")[,-1])
model1_sigma3_bayes2=as.matrix(read.csv("model1_sigma3_bayes2.csv")[,-1])

model2_sigma1_bayes1=as.matrix(read.csv("model2_sigma1_bayes1.csv")[,-1])
model2_sigma1_bayes2=as.matrix(read.csv("model2_sigma1_bayes2.csv")[,-1])
model2_sigma2_bayes1=as.matrix(read.csv("model2_sigma2_bayes1.csv")[,-1])
model2_sigma2_bayes2=as.matrix(read.csv("model2_sigma2_bayes2.csv")[,-1])
model2_sigma3_bayes1=as.matrix(read.csv("model2_sigma3_bayes1.csv")[,-1])
model2_sigma3_bayes2=as.matrix(read.csv("model2_sigma3_bayes2.csv")[,-1])

model3_sigma1_bayes1=as.matrix(read.csv("model3_sigma1_bayes1.csv")[,-1])
model3_sigma1_bayes2=as.matrix(read.csv("model3_sigma1_bayes2.csv")[,-1])
model3_sigma2_bayes1=as.matrix(read.csv("model3_sigma2_bayes1.csv")[,-1])
model3_sigma2_bayes2=as.matrix(read.csv("model3_sigma2_bayes2.csv")[,-1])
model3_sigma3_bayes1=as.matrix(read.csv("model3_sigma3_bayes1.csv")[,-1])
model3_sigma3_bayes2=as.matrix(read.csv("model3_sigma3_bayes2.csv")[,-1])

model4_sigma1_bayes1=as.matrix(read.csv("model4_sigma1_bayes1.csv")[,-1])
model4_sigma1_bayes2=as.matrix(read.csv("model4_sigma1_bayes2.csv")[,-1])
model4_sigma2_bayes1=as.matrix(read.csv("model4_sigma2_bayes1.csv")[,-1])
model4_sigma2_bayes2=as.matrix(read.csv("model4_sigma2_bayes2.csv")[,-1])
model4_sigma3_bayes1=as.matrix(read.csv("model4_sigma3_bayes1.csv")[,-1])
model4_sigma3_bayes2=as.matrix(read.csv("model4_sigma3_bayes2.csv")[,-1])

performance_gl=function(gamma,model){
  criterion=edge_count(solve(model))
  n=ncol(gamma)
  tp=numeric(length=n); tn=numeric(length=n); fp=numeric(length=n); fn=numeric(length=n)
  
  for(i in 1:n){
    temp=as.vector(gamma[,i])
    tp[i]=length(intersect(which(temp==1),which(criterion==1)))
    tn[i]=length(intersect(which(temp==0),which(criterion==0)))
    fp[i]=length(intersect(which(temp==1),which(criterion==0)))
    fn[i]=length(intersect(which(temp==0),which(criterion==1)))
  }
  
  sp=tn/(tn+fp)
  se=tp/(tp+fn)
  mcc=(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  sp=c(mean(sp),sd(sp))
  se=c(mean(se),sd(se))
  mcc=c(mean(mcc),sd(mcc))
  
  return(list(sp,se,mcc))
}

performance_gl(model1_sigma1_bayes1,model1_sigma1)
performance_gl(model1_sigma1_bayes2,model1_sigma1)
performance_gl(model1_sigma2_bayes1,model1_sigma2)
performance_gl(model1_sigma2_bayes2,model1_sigma2)
performance_gl(model1_sigma3_bayes1,model1_sigma3)
performance_gl(model1_sigma3_bayes2,model1_sigma3)

performance_gl(model2_sigma1_bayes1,model2_sigma1)
performance_gl(model2_sigma1_bayes2,model2_sigma1)
performance_gl(model2_sigma2_bayes1,model2_sigma2)
performance_gl(model2_sigma2_bayes2,model2_sigma2)
performance_gl(model2_sigma3_bayes1,model2_sigma3)
performance_gl(model2_sigma3_bayes2,model2_sigma3)

performance_gl(model3_sigma1_bayes1,model3_sigma1)
performance_gl(model3_sigma1_bayes2,model3_sigma1)
performance_gl(model3_sigma2_bayes1,model3_sigma2)
performance_gl(model3_sigma2_bayes2,model3_sigma2)
performance_gl(model3_sigma3_bayes1,model3_sigma3)
performance_gl(model3_sigma3_bayes2,model3_sigma3)

performance_gl(model4_sigma1_bayes1,model4_sigma1)
performance_gl(model4_sigma1_bayes2,model4_sigma1)
performance_gl(model4_sigma2_bayes1,model4_sigma2)
performance_gl(model4_sigma2_bayes2,model4_sigma2)
performance_gl(model4_sigma3_bayes1,model4_sigma3)
performance_gl(model4_sigma3_bayes2,model4_sigma3)
