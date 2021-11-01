source("performance.R")
source("model_choice&measure_ver2.R")
library(mvtnorm); library(spcov)

#AR1 model
#p=50, n=100(done)

ar1=diag(50)
for(i in 1:49){
  for(j in (i+1):50){
    ar1[i,j]=0.75^(abs(i-j))
  }
}
ar1=t(ar1)+ar1-diag(diag(ar1))
ar1=solve(ar1)

sssl_perf=bayesian_measure(100,1,100*0.06,0.1^3,12000,3000,4/49,1,ar1)

#MPP
sssl_perf[[1]];sssl_perf[[2]];sssl_perf[[3]];sssl_perf[[4]];sssl_perf[[5]];sssl_perf[[6]]
#0.1358505; 1.210301; 2.346749; 0.997449; 1; NA

#MAP
sssl_perf[[7]];sssl_perf[[8]];sssl_perf[[9]];sssl_perf[[10]];sssl_perf[[11]];sssl_perf[[12]]
#0.138858; 1.310162; 2.436117; 0.9863946; 1; NA

#p=100, n=200(done)

ar1=diag(100)
for(i in 1:99){
  for(j in (i+1):100){
    ar1[i,j]=0.75^(abs(i-j))
  }
}
ar1=t(ar1)+ar1-diag(diag(ar1))
ar1=solve(ar1)

sssl_perf=bayesian_measure(200,1,200*0.06,0.1^3,12000,3000,4/99,1,ar1)
#MPP
sssl_perf[[1]];sssl_perf[[2]];sssl_perf[[3]];sssl_perf[[4]];sssl_perf[[5]];sssl_perf[[6]]
#0.1515961; 1.714286; 3.544179; 0.9806226; 0.7979798; NA

#MAP
sssl_perf[[7]];sssl_perf[[8]];sssl_perf[[9]];sssl_perf[[10]];sssl_perf[[11]];
#0.1419682; 1.714286; 3.602841; 0.9800041; 0.8484848

#AR2 model
#p=50, n=100(done)
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
sssl_perf=bayesian_measure(100,1,100*0.06,0.1^3,12000,3000,4/49,1,ar2)
sssl_perf[[1]];sssl_perf[[2]];sssl_perf[[3]];sssl_perf[[4]];sssl_perf[[5]];sssl_perf[[6]]
#0.05718284; 0.5; 1.066903; 0.9867021 ;0.4948454

sssl_perf[[7]];sssl_perf[[8]];sssl_perf[[9]];sssl_perf[[10]];sssl_perf[[11]];sssl_perf[[12]]
#0.05718284; 0.5; 1.066903; 0.998227; 0.6907216

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
sssl_perf=bayesian_measure(200,1,200*0.06,0.1^3,10000,5000,4/99,1,ar2)
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