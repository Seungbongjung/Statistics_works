source("performance.R")
source("model_choice&measure.R")
library(mvtnorm)

#scenario2 (p=100, n=200)

random1=read.csv("random_structure1_100.csv")
random1=as.matrix(random1[,-1])
true_edge=edge_ind(random1,0.1^3)

rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

set.seed(2)
for(i in 1:100){
  samp=rmvnorm(200,sigma=random1)
  s=t(samp)%*%samp/200
  s_edge=edge_ind(s,0.1^3)
  
  rmse_samp[i]=rmse(s,random1)
  mnorm_samp[i]=mnorm(s,random1)
  spec_samp[i]=spec_norm(s,random1)
  
  acc=accuracy(s_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.07323958; 0.001793808

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.6728289; 0.1389011

#2norm
mean(spec_samp);sd(spec_samp)
#2.356923; 0.2203503

#specificity
mean(sp_samp); sd(sp_samp)
#0.02428357; 0.002141966

#sensitivity
mean(se_samp); sd(se_samp)
#0.9822818; 0.004107695

#mcc
mean(mcc_samp,na.rm=TRUE); sd(mcc_samp,na.rm=TRUE)
#NaN; NA

#scenario3 (p=50, n=100)

random2=read.csv("random_structure2_50.csv")
random2=as.matrix(random2[,-1])
true_edge=edge_ind(random2,0.1^3)

rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

set.seed(2)
for(i in 1:100){
  samp=rmvnorm(100,sigma=random2)
  s=t(samp)%*%samp/100
  s_edge=edge_ind(s,0.1^3)
  
  rmse_samp[i]=rmse(s,random2)
  mnorm_samp[i]=mnorm(s,random2)
  spec_samp[i]=spec_norm(s,random2)
  
  acc=accuracy(s_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.2563302; 0.007666335

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.9792487; 0.1448334

#2norm
mean(spec_samp);sd(spec_samp)
#4.723616; 0.4555399

#specificity
mean(sp_samp); sd(sp_samp)
#0.003208333; 0.001564077

#sensitivity
mean(se_samp); sd(se_samp)
#1;0

#mcc
mean(mcc_samp,na.rm=TRUE); sd(mcc_samp,na.rm=TRUE)
#NaN; NA

#scenario4 (p=50, n=100)

random3=read.csv("random_structure3_50.csv")
random3=as.matrix(random3[,-1])
true_edge=edge_ind(random3,0.1^3)

rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

set.seed(2)
for(i in 1:100){
  samp=rmvnorm(100,sigma=random3)
  s=t(samp)%*%samp/100
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
#0.08422206; 0.003008909

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.3136207; 0.04283796

#2norm
mean(spec_samp);sd(spec_samp)
#1.676928; 0.1785506

#specificity
mean(sp_samp); sd(sp_samp)
#0.01010204; 0.002578869

#sensitivity
mean(se_samp); sd(se_samp)
#1;0

#mcc
mean(mcc_samp,na.rm=TRUE); sd(mcc_samp,na.rm=TRUE)
#NaN; NA

#scenario5 (p=100, n=200)

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
  samp=rmvnorm(200,sigma=model)
  s=t(samp)%*%samp/200
  s_edge=edge_ind(s,0.1^3)
  
  rmse_samp[i]=rmse(s,model)
  mnorm_samp[i]=mnorm(s,model)
  spec_samp[i]=spec_norm(s,model)
  
  acc=accuracy(s_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.07118934; 0.001209381

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.2970257; 0.03589285

#2norm
mean(spec_samp);sd(spec_samp)
#2.078636; 0.1585887

#specificity
mean(sp_samp); sd(sp_samp)
#0.01103692; 0.001693901

#sensitivity
mean(se_samp); sd(se_samp)
#0.9977711; 0.001416878

#mcc
mean(mcc_samp,na.rm=TRUE); sd(mcc_samp,na.rm=TRUE)
#NaN; NA


#scenario6 (p=50,n=100)

toep=matrix(rep(0,50^2),ncol=50)
for(i in 1:49){
  for(j in (i+1):50){
    toep[i,j]=0.75^(abs(i-j))
  }
}
toep=t(toep)+toep+diag(rep(1,50))

toep=solve(toep)
true_edge=edge_ind(toep,0.1^3)
sum(true_edge)
#sample covariance
rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

for(i in 1:100){
  samp=rmvnorm(100,sigma=toep)
  s=t(samp)%*%samp/100
  s_edge=edge_ind(s,0.1^3)
  
  rmse_samp[i]=rmse(s,toep)
  mnorm_samp[i]=mnorm(s,toep)
  spec_samp[i]=spec_norm(s,toep)
  
  acc=accuracy(s_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.356205; 0.01385241

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#1.337605; 0.1715263

#2norm
mean(spec_samp);sd(spec_samp)
#7.128174; 0.6934294

#specificity
mean(sp_samp); sd(sp_samp)
# 0.002363946; 0.001375812

#sensitivity
mean(se_samp); sd(se_samp)
#1; 0

#mcc
mean(mcc_samp,na.rm=TRUE); sd(mcc_samp,na.rm=TRUE)
#NaN; NA
