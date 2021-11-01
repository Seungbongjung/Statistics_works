source("model_choice&measure.R")
library(glasso)

#scenario5
toep=matrix(rep(0,100^2),ncol=100)
for(i in 1:99){
  for(j in (i+1):100){
    toep[i,j]=0.75^(abs(i-j))
  }
}
toep=t(toep)+toep+diag(rep(1,100))
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
#0.2140305; 0.006353871

#mnorm
mean(mnorm_samp1); sd(mnorm_samp1)
#0.8902245; 0.1138324

#2norm
mean(spec_samp1);sd(spec_samp1)
#7.549499; 0.8459087

#specificity
mean(sp_samp1); sd(sp_samp1)
#  0.002852792; 0.0017055

#sensitivity
mean(se_samp1); sd(se_samp1)
#0.003783784;  0.001117014

#mcc
mean(mcc_samp1,na.rm=TRUE); sd(mcc_samp1,na.rm=TRUE)
#NaN; NA

#scenario6
banded=matrix(rep(0,100^2),ncol=100)
for(i in 1:100){
  for(j in i:100){
    if(abs(i-j)<=20){
      banded[i,j]=1-abs(i-j)/20
    }
  }
}
banded=t(banded)+banded-diag(diag(banded))
true_edge=edge_ind(banded,0.1^3)

#sample covariance
rmse_samp2=numeric(length=100); mnorm_samp2=numeric(length=100); spec_samp2=numeric(length=100)
sp_samp2=numeric(length=100); se_samp2=numeric(length=100); mcc_samp2=numeric(length=100)

for(i in 1:100){
  samp=rmvnorm(200,sigma=banded)
  s=t(samp)%*%samp/200
  s_edge=edge_ind(s,0.1^3)
  
  rmse_samp2[i]=rmse(s,banded)
  mnorm_samp2[i]=mnorm(s,banded)
  spec_samp2[i]=spec_norm(s,banded)
  
  acc=accuracy(s_edge,true_edge)
  sp_samp2[i]=acc[1]
  se_samp2[i]=acc[2]
  mcc_samp2[i]=acc[3]
}

#rmse
mean(rmse_samp2); sd(rmse_samp2)
#0.07429439; 0.01015552

#mnorm
mean(mnorm_samp2); sd(mnorm_samp2)
#0.2553541; 0.04063939

#2norm
mean(spec_samp2);sd(spec_samp2)
#4.752106; 1.048046

#specificity
mean(sp_samp2); sd(sp_samp2)
# 0.01130247; 0.002807849

#sensitivity
mean(se_samp2); sd(se_samp2)
#0.9993684; 0.0006938337

#mcc
mean(mcc_samp2,na.rm=TRUE); sd(mcc_samp2,na.rm=TRUE)
#NaN; NA
