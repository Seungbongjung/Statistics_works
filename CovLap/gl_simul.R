source("performance.R")
source("model_choice&measure_ver2.R")
library(mvtnorm)
library(spcov)

#scenario2 (p=100, n=200)

step.size=100
tol=0.001
P=matrix(1,100,100)
diag(P)=0
lam=0.06

random1=read.csv("random_structure1_100.csv")
random1=as.matrix(random1[,-1])
true_edge=edge_ind(random1,0.1^3)

rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

set.seed(2)
for(i in 1:100){
  print(paste("replication: ",i,sep=""))
  
  samp=rmvnorm(200,sigma=random1)
  s=t(samp)%*%samp/200
  
  mm=spcov(Sigma=s,S=s,lambda=lam*P,step.size=step.size,n.inner.steps=200,
           thr.inner=0,tol.outer=tol,trace=1)$Sigma
  
  gl_edge=edge_ind(mm,0.1^3)
  
  rmse_samp[i]=rmse(mm,random1)
  mnorm_samp[i]=mnorm(mm,random1)
  spec_samp[i]=spec_norm(mm,random1)
  
  acc=accuracy(gl_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.04097324; 0.001668946

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.9256504; 0.1828577

#2norm
mean(spec_samp);sd(spec_samp)
#1.263337; 0.1307285

#specificity
mean(sp_samp); sd(sp_samp)
#0.3486522; 0.009242691

#sensitivity
mean(se_samp); sd(se_samp)
#0.6545846; 0.0151022

#mcc
mean(mcc_samp,na.rm=TRUE); sd(mcc_samp,na.rm=TRUE)
#NaN; NA

#scenario3 (p=50, n=100)

step.size=100
tol=0.001
P=matrix(1,50,50)
diag(P)=0
lam=0.06

random2=read.csv("random_structure2_50.csv")
random2=as.matrix(random2[,-1])
true_edge=edge_ind(random2,0.1^3)

rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

set.seed(2)
for(i in 1:100){
  print(paste(c("replication: ",i),sep=""))
  
  samp=rmvnorm(100,sigma=random2)
  s=t(samp)%*%samp/100
  
  mm=spcov(Sigma=s,S=s,lambda=lam*P,step.size=step.size,n.inner.steps=200,
           thr.inner=0,tol.outer=tol,trace=1)$Sigma
  
  gl_edge=edge_ind(mm,0.1^3)
  
  rmse_samp[i]=rmse(mm,random2)
  mnorm_samp[i]=mnorm(mm,random2)
  spec_samp[i]=spec_norm(mm,random2)
  
  acc=accuracy(gl_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.1390417; 0.006674287

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#1.171516; 0.09401429

#2norm
mean(spec_samp);sd(spec_samp)
#2.434382; 0.1700338

#specificity
mean(sp_samp); sd(sp_samp)
#0.5833833; 0.02651996

#sensitivity
mean(se_samp); sd(se_samp)
#0.9996;  0.004

#mcc
mean(mcc_samp,na.rm=TRUE); sd(mcc_samp,na.rm=TRUE)
#NaN; NA

#scenario4 (p=50, n=100)

step.size=100
tol=0.001
P=matrix(1,50,50)
diag(P)=0
lam=0.06

random3=read.csv("random_structure3_50.csv")
random3=as.matrix(random3[,-1])
true_edge=edge_ind(random3,0.1^3)

rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

set.seed(2)
for(i in 1:100){
  print(paste(c("replication: ",i),sep=""))
  
  samp=rmvnorm(100,sigma=random3)
  s=t(samp)%*%samp/100
  
  mm=spcov(Sigma=s,S=s,lambda=lam*P,step.size=step.size,n.inner.steps=200,
           thr.inner=0,tol.outer=tol,trace=1)$Sigma
  
  gl_edge=edge_ind(mm,0.1^3)
  
  rmse_samp[i]=rmse(mm,random3)
  mnorm_samp[i]=mnorm(mm,random3)
  spec_samp[i]=spec_norm(mm,random3)
  
  acc=accuracy(gl_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.05371842; 0.001841355

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.3512786; 0.03439454

#2norm
mean(spec_samp);sd(spec_samp)
#0.9658245; 0.05545525

#specificity
mean(sp_samp); sd(sp_samp)
#0.235102; 0.01377217

#sensitivity
mean(se_samp); sd(se_samp)
#1;0

#mcc
mean(mcc_samp,na.rm=TRUE); sd(mcc_samp,na.rm=TRUE)
#NaN; NA

#scenario5 (p=100, n=200)

step.size=100
tol=0.01
P=matrix(1,100,100)
diag(P)=0
lam=0.06

model=matrix(rep(0,100^2),ncol=100)
for(i in 1:99){
  for(j in (i+1):100){
      model[i,j]=sample(c(0.5,0),size=1,prob=c(0.2,0.8))
  }
}
model=model+t(model); diag(model)=1
delta=max(-min(eigen(model)$values),0)+0.05
model=(model+delta*diag(100))/(1+delta)
true_edge=edge_ind(model,0.1^3)

#sample covariance
rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

for(i in 1:100){
  print(paste(c("replication: ",i),sep=""))
  
  samp=rmvnorm(200,sigma=model)
  s=t(samp)%*%samp/200
  
  mm=spcov(Sigma=s,S=s,lambda=lam*P,step.size=step.size,n.inner.steps=200,
           thr.inner=0,tol.outer=tol,trace=1)$Sigma
  
  gl_edge=edge_ind(mm,0.1^3)
  
  rmse_samp[i]=rmse(mm,model)
  mnorm_samp[i]=mnorm(mm,model)
  spec_samp[i]=spec_norm(mm,model)
  
  acc=accuracy(gl_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.05092924; 0.0007512975

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.422388; 0.02450827

#2norm
mean(spec_samp);sd(spec_samp)
#1.974601; 0.07835921

#specificity
mean(sp_samp); sd(sp_samp)
#0.3392135; 0.009460125

#sensitivity
mean(se_samp); sd(se_samp)
#0.9437952; 0.007301654

#mcc
mean(mcc_samp,na.rm=TRUE); sd(mcc_samp,na.rm=TRUE)
#NaN; NA


#scenario6 (p=50,n=100)

step.size=100
tol=0.001
P=matrix(1,50,50)
diag(P)=0
lam=0.06

toep=matrix(rep(0,50^2),ncol=50)
for(i in 1:49){
  for(j in (i+1):50){
    toep[i,j]=0.75^(abs(i-j))
  }
}
toep=t(toep)+toep+diag(rep(1,50))
toep=solve(toep)
true_edge=edge_ind(toep,0.1^3)

#sample covariance
rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

for(i in 1:100){
  print(paste(c("replication: ",i),sep=""))
  
  samp=rmvnorm(100,sigma=toep)
  s=t(samp)%*%samp/100
  
  mm=spcov(Sigma=s,S=s,lambda=lam*P,step.size=step.size,n.inner.steps=200,
           thr.inner=0,tol.outer=tol,trace=1)$Sigma
  
  gl_edge=edge_ind(mm,0.1^3)
  
  rmse_samp[i]=rmse(mm,toep)
  mnorm_samp[i]=mnorm(mm,toep)
  spec_samp[i]=spec_norm(mm,toep)
  
  acc=accuracy(gl_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.2101934; 0.009052997

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#1.705604; 0.1492609

#2norm
mean(spec_samp);sd(spec_samp)
#3.334658; 0.2316332

#specificity
mean(sp_samp); sd(sp_samp)
# 0.5886905; 0.01999217

#sensitivity
mean(se_samp); sd(se_samp)
#1; 0

#mcc
mean(mcc_samp,na.rm=TRUE); sd(mcc_samp,na.rm=TRUE)
#NaN; NA
