source("model_choice&measure.R")
library(mvtnorm)

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

set.seed(1)
for(i in 1:100){
  samp=rmvnorm(250,sigma=curr_ex)
  s=t(samp)%*%samp/250
  s_edge=edge_ind(s,0.1^3)
  
  rmse_samp[i]=rmse(s,curr_ex)
  mnorm_samp[i]=mnorm(s,curr_ex)
  spec_samp[i]=spec_norm(s,curr_ex)
  
  acc=accuracy(s_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}
#rmse
mean(rmse_samp); sd(rmse_samp)
#0.0297044; 0.005607053

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.1539036; 0.07603278

#2norm
mean(spec_samp);sd(spec_samp)
#0.2528385; 0.07936456

#specificity
mean(sp_samp); sd(sp_samp)
# 0.04188679; 0.02957009

#sensitivity
mean(se_samp); sd(se_samp)
#0.9838462; 0.03508062

#mcc
mean(mcc_samp,na.rm=TRUE); sd(mcc_samp,na.rm=TRUE)
#0.05591463; 0.09171744

#scenario2-random structure1
random1=read.csv("random_structure1.csv")
random1=as.matrix(random1[,-1])
true_edge=edge_ind(random1,0.1^3)

rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

set.seed(2)
for(i in 1:100){
  samp=rmvnorm(100,sigma=random1)
  s=t(samp)%*%samp/100
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
#0.1101033; 0.00468437

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.8125537; 0.1915116

#2norm
mean(spec_samp);sd(spec_samp)
#2.31521; 0.2776134

#specificity
mean(sp_samp); sd(sp_samp)
# 0.01679188; 0.004334765

#sensitivity
mean(se_samp); sd(se_samp)
#0.9894583;  0.006308091

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
for(i in 1:100){
  samp=rmvnorm(200,sigma=random2)
  s=t(samp)%*%samp/200
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
#0.2221888; 0.003770648

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.9156252; 0.1012726

#2norm
mean(spec_samp);sd(spec_samp)
#6.128957; 0.3786077

#specificity
mean(sp_samp); sd(sp_samp)
# 0.003540552;  0.0008728405

#sensitivity
mean(se_samp); sd(se_samp)
#1;  0

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
#0.05838504; 0.001177863

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.2374097; 0.2374097

#2norm
mean(spec_samp);sd(spec_samp)
#1.73451; 0.1002527

#specificity
mean(sp_samp); sd(sp_samp)
# 0.01395589; 0.001611095

#sensitivity
mean(se_samp); sd(se_samp)
#1;  0

#mcc
mean(mcc_samp,na.rm=TRUE); sd(mcc_samp,na.rm=TRUE)
#NaN; NA