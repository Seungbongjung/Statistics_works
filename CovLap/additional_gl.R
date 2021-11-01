source("performance.R")
library(mvtnorm); library(spcov)

#AR1 model

step.size=100
tol=0.001
P=matrix(1,50,50)
diag(P)=0
lam=0.06

#p=50, n=100

ar1=diag(50)
for(i in 1:49){
  for(j in (i+1):50){
    ar1[i,j]=0.75^(abs(i-j))
  }
}
ar1=t(ar1)+ar1-diag(diag(ar1))
ar1=solve(ar1)
true_edge=edge_ind(ar1,0.1^3)

rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

for(i in 1:100){
  samp=rmvnorm(n=100,sigma=ar1)
  s=t(samp)%*%samp/100
  
  print(paste(i,"th replication",sep=""))
  mm=spcov(Sigma=s,S=s,lambda=lam*P,step.size=step.size,n.inner.steps=200,
           thr.inner=0,tol.outer=tol,trace=1)$Sigma
  
  gl_edge=edge_ind(mm,0.1^3)
  
  rmse_samp[i]=rmse(mm,ar1)
  mnorm_samp[i]=mnorm(mm,ar1)
  spec_samp[i]=spec_norm(mm,ar1)
  
  acc=accuracy(gl_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.2107166; 0.009440502

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#1.6986; 0.1146048

#2norm
mean(spec_samp); sd(spec_samp)
#3.339946; 0.2166977

#specificity 
mean(sp_samp); sd(sp_samp)
#0.5889881; 0.0171596

#sensitivity
mean(se_samp); sd(se_samp)
# 1;0

#mcc
mean(mcc_samp); sd(mcc_samp)
#NaN; NA

#p=100, n=200

step.size=100
tol=0.001
P=matrix(1,100,100)
diag(P)=0
lam=0.06

ar1=diag(100)
for(i in 1:99){
  for(j in (i+1):100){
    ar1[i,j]=0.75^(abs(i-j))
  }
}
ar1=t(ar1)+ar1-diag(diag(ar1))
ar1=solve(ar1)
true_edge=edge_ind(ar1,0.1^3)

rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

for(i in 1:100){
  samp=rmvnorm(n=200,sigma=ar1)
  s=t(samp)%*%samp/200
  
  print(paste(i,"th replication",sep=""))
  mm=spcov(Sigma=s,S=s,lambda=lam*P,step.size=step.size,n.inner.steps=200,
           thr.inner=0,tol.outer=tol,trace=1)$Sigma
  
  gl_edge=edge_ind(mm,0.1^3)
  
  rmse_samp[i]=rmse(mm,ar1)
  mnorm_samp[i]=mnorm(mm,ar1)
  spec_samp[i]=spec_norm(mm,ar1)
  
  acc=accuracy(gl_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.1291366;  0.003564975

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#1.498024; 0.09159307

#2norm
mean(spec_samp); sd(spec_samp)
#2.828411; 0.12583

#specificity 
mean(sp_samp); sd(sp_samp)
#0.7258957; 0.008913763

#sensitivity
mean(se_samp); sd(se_samp)
# 1;0

#mcc
mean(mcc_samp); sd(mcc_samp)
#NA; NA

#AR2 model
#p=50, n=100
step.size=100
tol=0.001
P=matrix(1,50,50)
diag(P)=0
lam=0.06

ar2=diag(50)
for(i in 1:49){
  for(j in (i+1):50){
    if((j-i)==1){
      ar2[i,j]=0.5
    }
    if((j-i)==2){
      ar2[i,j]=0.25
    }
  }
}
ar2=t(ar2)+ar2-diag(diag(ar2))
true_edge=edge_ind(ar2,0.1^3)

rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

for(i in 1:100){
  samp=rmvnorm(n=100,sigma=ar2)
  s=t(samp)%*%samp/100
  
  print(paste(i,"th replication",sep=""))
  mm=spcov(Sigma=s,S=s,lambda=lam*P,step.size=step.size,n.inner.steps=200,
           thr.inner=0,tol.outer=tol,trace=1)$Sigma
  
  gl_edge=edge_ind(mm,0.1^3)
  
  rmse_samp[i]=rmse(mm,ar2)
  mnorm_samp[i]=mnorm(mm,ar2)
  spec_samp[i]=spec_norm(mm,ar2)
  
  acc=accuracy(gl_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.07305428; 0.002517245

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.4473803; 0.03652329

#2norm
mean(spec_samp); sd(spec_samp)
#1.353454; 0.07890297

#specificity 
mean(sp_samp); sd(sp_samp)
#0.2484929; 0.01229019

#sensitivity
mean(se_samp); sd(se_samp)
# 0.9945361;  0.007816367

#mcc
mean(mcc_samp,na.rm=TRUE); sd(mcc_samp,na.rm=TRUE)
#NA; NA

#p=100, n=200
step.size=100
tol=0.001
P=matrix(1,100,100)
diag(P)=0
lam=0.06


ar2=diag(100)
for(i in 1:99){
  for(j in (i+1):100){
    if((j-i)==1){
      ar2[i,j]=0.5
    }
    if((j-i)==2){
      ar2[i,j]=0.25
    }
  }
}
ar2=t(ar2)+ar2-diag(diag(ar2))
true_edge=edge_ind(ar2,0.1^3)

rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

for(i in 1:100){
  samp=rmvnorm(n=200,sigma=ar2)
  s=t(samp)%*%samp/200
  
  print(paste(i,"th replication",sep=""))
  mm=spcov(Sigma=s,S=s,lambda=lam*P,step.size=step.size,n.inner.steps=200,
           thr.inner=0,tol.outer=tol,trace=1)$Sigma
  
  gl_edge=edge_ind(mm,0.1^3)
  
  rmse_samp[i]=rmse(mm,ar2)
  mnorm_samp[i]=mnorm(mm,ar2)
  spec_samp[i]=spec_norm(mm,ar2)
  
  acc=accuracy(gl_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.05004387; 0.001201792

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.4237814; 0.02236823

#2norm
mean(spec_samp); sd(spec_samp)
#1.309292; 0.04836717

#specificity 
mean(sp_samp); sd(sp_samp)
#0.3258679; 0.007975205

#sensitivity
mean(se_samp); sd(se_samp)
# 0.9996447; 0.001301687

#mcc
mean(mcc_samp); sd(mcc_samp)
#NA; NA



