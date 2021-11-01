source("model_choice&measure.R")
library(mvtnorm)
library(spcov)

#scenario1
#sparse covariance that mimics currency exchange rate, n=250
curr_ex=matrix(rep(0,12^2),nrow=12)
diag(curr_ex)=c(0.239,1.554,0.362,0.199,0.349,0.295,0.715,0.164,0.518,0.379,0.159,0.207)
curr_ex[1,2]=0.117; curr_ex[1,8]=0.031; curr_ex[3,4]=0.002; curr_ex[4,5]=0.094
curr_ex[5,12]=-0.036; curr_ex[6,7]=-0.229; curr_ex[6,8]=0.002; curr_ex[8,9]=0.112
curr_ex[8,10]=-0.028; curr_ex[8,11]=-0.008; curr_ex[9,10]=-0.193; curr_ex[9,11]=-0.090
curr_ex[10,11]=0.167

curr_ex=curr_ex+t(curr_ex)-diag(diag(curr_ex))
true_edge=edge_ind(curr_ex,0.1^3)

rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

step.size=100
tol=1e-3
P=matrix(1,12,12)
diag(P)=0
lam=0.06

for(i in 1:100){
  samp=rmvnorm(n=250,sigma=curr_ex)
  s=t(samp)%*%samp/250
  
  mm=spcov(Sigma=s,S=s,lambda=lam*P,step.size=step.size,n.inner.steps=200,
           thr.inner=0,tol.outer=tol,trace=1)$Sigma
  
  gl_edge=edge_ind(mm,0.1^3)
  
  rmse_samp[i]=rmse(s,curr_ex)
  mnorm_samp[i]=mnorm(s,curr_ex)
  spec_samp[i]=spec_norm(s,curr_ex)
  
  acc=accuracy(gl_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}
#rmse
mean(rmse_samp); sd(rmse_samp)
#0.02882296; 0.004588071

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.1381143; 0.06853715

#2norm
mean(spec_samp);sd(spec_samp)
#0.2369619; 0.06330581

#specificity
mean(sp_samp); sd(sp_samp)
# 0.2883019; 0.06431106

#sensitivity
mean(se_samp); sd(se_samp)
#0.9615385; 0.05069595

#mcc
mean(mcc_samp); sd(mcc_samp)
#0.2320372; 0.07412892

#scenario2-random structure1
random1=read.csv("random_structure1.csv")
random1=as.matrix(random1[,-1])
true_edge=edge_ind(random1,0.1^3)

rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

step.size=100
tol=1e-3
P=matrix(1,50,50)
diag(P)=0
lam=0.06

set.seed(2)
for(i in 1:100){
  samp=rmvnorm(100,sigma=random1)
  s=t(samp)%*%samp/100
  
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
#0.07753903; 0.004016276

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.7973847; 0.1683615

#2norm
mean(spec_samp);sd(spec_samp)
#1.482957; 0.1709074

#specificity
mean(sp_samp); sd(sp_samp)
#0.2519695; 0.01743604

#sensitivity
mean(se_samp); sd(se_samp)
#0.77; 0.02737716

#mcc
mean(mcc_samp,na.rm=TRUE); sd(mcc_samp,na.rm=TRUE)
#NaN; NA

#scenario3-random structure2
random2=read.csv("random_structure2.csv",header=FALSE)
random2=as.matrix(random2)
true_edge=edge_ind(random2,0.1^3)

rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

set.seed(3)

step.size=100
tol=1e-3
P=matrix(1,100,100)
diag(P)=0
lam=0.06

for(i in 1:100){
  samp=rmvnorm(200,sigma=random2)
  s=t(samp)%*%samp/200
  
  print(paste(c("iteration: ",i),sep=""))
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
#0.08437705; 0.003049868

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#1.169109;  0.08552816

#2norm
mean(spec_samp);sd(spec_samp)
#2.249376; 0.1474289

#specificity
mean(sp_samp); sd(sp_samp)
# 0.899576;  0.008369846

#sensitivity
mean(se_samp); sd(se_samp)
#0.9995652; 0.002140722

#mcc
mean(mcc_samp,na.rm=TRUE); sd(mcc_samp,na.rm=TRUE)
#NaN; NA

#scenario4-random structure3
random3=read.csv("random_structure3.csv",header=FALSE)
random3=as.matrix(random3)
true_edge=edge_ind(random3,0.1^3)

rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

set.seed(3)
for(i in 1:100){
  samp=rmvnorm(200,sigma=random3)
  s=t(samp)%*%samp/200
  s_edge=edge_ind(s,0.1^3)
  
  rmse_samp[i]=rmse(s,random3)
  mnorm_samp[i]=mnorm(s,random3)
  spec_samp[i]=spec_norm(s,random3)
  
  acc=accuracy(s_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.03404026; 0.0007835315

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.3192972; 0.02367287

#2norm
mean(spec_samp);sd(spec_samp)
#0.90000059; 0.03482349

#specificity
mean(sp_samp); sd(sp_samp)
# 0.6579365; 0.007261079

#sensitivity
mean(se_samp); sd(se_samp)
#1;  0

#mcc
mean(mcc_samp,na.rm=TRUE); sd(mcc_samp,na.rm=TRUE)
#NaN; NA

