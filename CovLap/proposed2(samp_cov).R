source("model_choice&measure.R")

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

#sample covariance
rmse_samp1=numeric(length=100); mnorm_samp1=numeric(length=100); spec_samp1=numeric(length=100)
sp_samp1=numeric(length=100); se_samp1=numeric(length=100); mcc_samp1=numeric(length=100)

for(i in 1:100){
  samp=rmvnorm(200,sigma=toep)
  s=t(samp)%*%samp/200
  s_edge=edge_ind(s,0.1^3)

  rmse_samp1[i]=rmse(s,toep)
  mnorm_samp1[i]=mnorm(s,toep)
  spec_samp1[i]=spec_norm(s,toep)
  
  acc=accuracy(s_edge,true_edge)
  sp_samp1[i]=acc[1]
  se_samp1[i]=acc[2]
  mcc_samp1[i]=acc[3]
}

#rmse
mean(rmse_samp1); sd(rmse_samp1)
#0.2526401; 0.004469016

#mnorm
mean(mnorm_samp1); sd(mnorm_samp1)
#1.046026; 0.09237194

#2norm
mean(spec_samp1);sd(spec_samp1)
#7.554739; 0.4896099

#specificity
mean(sp_samp1); sd(sp_samp1)
# 0.003160173; 0.0008172808

#sensitivity
mean(se_samp1); sd(se_samp1)
#1; 0

#mcc
mean(mcc_samp1,na.rm=TRUE); sd(mcc_samp1,na.rm=TRUE)
#NaN; NA

#scenario6
true_edge=edge_ind(model,0.1^3)

#sample covariance
rmse_samp2=numeric(length=100); mnorm_samp2=numeric(length=100); spec_samp2=numeric(length=100)
sp_samp2=numeric(length=100); se_samp2=numeric(length=100); mcc_samp2=numeric(length=100)

for(i in 1:100){
  samp=rmvnorm(100,sigma=model)
  s=t(samp)%*%samp/100
  s_edge=edge_ind(s,0.1^3)
  
  rmse_samp2[i]=rmse(s,model)
  mnorm_samp2[i]=mnorm(s,model)
  spec_samp2[i]=spec_norm(s,model)
  
  acc=accuracy(s_edge,true_edge)
  sp_samp2[i]=acc[1]
  se_samp2[i]=acc[2]
  mcc_samp2[i]=acc[3]
}

#rmse
mean(rmse_samp2); sd(rmse_samp2)
#0.101335; 0.003495427

#mnorm
mean(mnorm_samp2); sd(mnorm_samp2)
#0.3857794; 0.05454149

#2norm
mean(spec_samp2);sd(spec_samp2)
#1.993569; 0.210729

#specificity
mean(sp_samp2); sd(sp_samp2)
# 0.008220676; 0.002969883

#sensitivity
mean(se_samp2); sd(se_samp2)
#0.9990411; 0.001978691

#mcc
mean(mcc_samp2,na.rm=TRUE); sd(mcc_samp2,na.rm=TRUE)
#0.02957149; 0.0153406

