library(spcov)
source("model_choice&measure.R")

step.size=100
tol=0.001
P=matrix(1,100,100)
diag(P)=0
lam=0.06

#scenario5
toep=matrix(rep(0,100^2),ncol=100)
for(i in 1:99){
  for(j in (i+1):100){
    toep[i,j]=0.75^(abs(i-j))
  }
}
toep=t(toep)+toep+diag(rep(1,100))
toep=solve(toep)
true_edge=edge_ind(toep,0.1^3)

#graphical lasso
rmse_samp1=numeric(length=100); mnorm_samp1=numeric(length=100); spec_samp1=numeric(length=100)
sp_samp1=numeric(length=100); se_samp1=numeric(length=100); mcc_samp1=numeric(length=100)

for(i in 1:100){
  print(paste(c("replication: ",i),sep=""))

  samp=rmvnorm(200,sigma=toep)
  s=t(samp)%*%samp/200
  
  mm=spcov(Sigma=s,S=s,lambda=lam*P,step.size=step.size,n.inner.steps=200,
           thr.inner=0,tol.outer=tol,trace=1)$Sigma
  
  gl_edge=edge_ind(mm,0.1^3)
  
  rmse_samp1[i]=rmse(mm,toep)
  mnorm_samp1[i]=mnorm(mm,toep)
  spec_samp1[i]=spec_norm(mm,toep)
  
  acc=accuracy(gl_edge,true_edge)
  sp_samp1[i]=acc[1]
  se_samp1[i]=acc[2]
  mcc_samp1[i]=acc[3]
}

#rmse
mean(rmse_samp1); sd(rmse_samp1)
#0.1064871; 0.001516619

#mnorm
mean(mnorm_samp1); sd(mnorm_samp1)
#0.5072663; 0.01928159

#2norm
mean(spec_samp1);sd(spec_samp1)
#4.576636; 0.0692268

#specificity
mean(sp_samp1); sd(sp_samp1)
# 0.2680772;  0.008535459

#sensitivity
mean(se_samp1); sd(se_samp1)
#0.7906905; 0.008151595

#mcc
mean(mcc_samp1,na.rm=TRUE); sd(mcc_samp1,na.rm=TRUE)
#NaN; NA

#scenario6
step.size=100
tol=0.001
P=matrix(1,50,50)
diag(P)=0
lam=0.06

model=matrix(rep(0,50^2),ncol=50)
s0=50
for(i in 1:49){
  for(j in (i+1):50){
    if(i<=s0){
      model[i,j]=sample(c(0.5,0),size=1,prob=c(0.2,0.8))
    }else{
      model[i,j]=0.5
    }
  }
}
model=model+t(model); diag(model)=1
delta=max(-min(eigen(model)$values),0)+0.05
model=(model+delta*diag(50))/(1+delta)

#graphical lasso
rmse_samp2=numeric(length=100); mnorm_samp2=numeric(length=100); spec_samp2=numeric(length=100)
sp_samp2=numeric(length=100); se_samp2=numeric(length=100); mcc_samp2=numeric(length=100)

for(i in 1:100){
  print(paste(c("replication: ",i),sep=""))
  
  samp=rmvnorm(100,sigma=model)
  s=t(samp)%*%samp/100
  
  mm=spcov(Sigma=s,S=s,lambda=lam*P,step.size=step.size,n.inner.steps=200,
           thr.inner=0,tol.outer=tol,trace=1)$Sigma
  
  gl_edge=edge_ind(mm,0.1^3)
  
  rmse_samp2[i]=rmse(mm,model)
  mnorm_samp2[i]=mnorm(mm,model)
  spec_samp2[i]=spec_norm(mm,model)
  
  acc=accuracy(gl_edge,true_edge)
  sp_samp2[i]=acc[1]
  se_samp2[i]=acc[2]
  mcc_samp2[i]=acc[3]
}

#rmse
mean(rmse_samp2); sd(rmse_samp2)
#0.07164835; 0.00193795

#mnorm
mean(mnorm_samp2); sd(mnorm_samp2)
#0.4426029; 0.03410894

#2norm
mean(spec_samp2);sd(spec_samp2)
#1.43289; 0.1078849

#specificity
mean(sp_samp2); sd(sp_samp2)
# 0.2622266; 0.0132523

#sensitivity
mean(se_samp2); sd(se_samp2)
#0.973379; 0.01062589

#mcc
mean(mcc_samp2,na.rm=TRUE); sd(mcc_samp2,na.rm=TRUE)
#NaN; NA
