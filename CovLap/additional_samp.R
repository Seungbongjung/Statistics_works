source("performance.R")
library(mvtnorm); library(spcov)

#AR1 model

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
  
  samp_edge=edge_ind(s,0.1^3)
  
  rmse_samp[i]=rmse(s,ar1)
  mnorm_samp[i]=mnorm(s,ar1)
  spec_samp[i]=spec_norm(s,ar1)
  
  acc=accuracy(samp_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.3553866; 0.01278761

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#1.365836; 0.206505

#2norm
mean(spec_samp); sd(spec_samp)
#7.155482; 0.8180765

#specificity 
mean(sp_samp); sd(sp_samp)
#0.002380952;  0.001261837

#sensitivity
mean(se_samp); sd(se_samp)
# 1;0

#mcc
mean(mcc_samp); sd(mcc_samp)
#NaN; NA

#p=100, n=200

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
  
  samp_edge=edge_ind(s,0.1^3)
  
  rmse_samp[i]=rmse(s,ar1)
  mnorm_samp[i]=mnorm(s,ar1)
  spec_samp[i]=spec_norm(s,ar1)
  
  acc=accuracy(samp_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.2515424; 0.004771637

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#1.061335; 0.1288696

#2norm
mean(spec_samp); sd(spec_samp)
#7.463462; 0.4440703

#specificity 
mean(sp_samp); sd(sp_samp)
#0.003302412; 0.0008968737

#sensitivity
mean(se_samp); sd(se_samp)
# 1;0

#mcc
mean(mcc_samp); sd(mcc_samp)
#NA; NA

#AR2 model
#p=50, n=100

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
  
  samp_edge=edge_ind(s,0.1^3)
  
  rmse_samp[i]=rmse(s,ar2)
  mnorm_samp[i]=mnorm(s,ar2)
  spec_samp[i]=spec_norm(s,ar2)
  
  acc=accuracy(samp_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.1019923; 0.00472153

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.3909498; 0.06123361

#2norm
mean(spec_samp); sd(spec_samp)
#2.144687; 0.2370766

#specificity 
mean(sp_samp); sd(sp_samp)
#0.008537234; 0.002512705

#sensitivity
mean(se_samp); sd(se_samp)
# 0.9997938; 0.00145057

#mcc
mean(mcc_samp); sd(mcc_samp)
#NA; NA

#p=100, n=200
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
  
  samp_edge=edge_ind(s,0.1^3)
  
  rmse_samp[i]=rmse(s,ar2)
  mnorm_samp[i]=mnorm(s,ar2)
  spec_samp[i]=spec_norm(s,ar2)
  
  acc=accuracy(samp_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.07098919; 0.001312084

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.291229; 0.0311699

#2norm
mean(spec_samp); sd(spec_samp)
#2.145154; 0.1345521

#specificity 
mean(sp_samp); sd(sp_samp)
#0.01126588; 0.001331638

#sensitivity
mean(se_samp); sd(se_samp)
# 1;0

#mcc
mean(mcc_samp); sd(mcc_samp)
#NA; NA

#circle model
#p=50, n=100
circle=2*diag(50)
for(i in 1:49){
  for(j in (i+1):50){
    if((j-i)==1){
      circle[i,j]=1
    }
    if(j==50){
      circle[1,j]=0.9
    }
  }
}
circle=t(circle)+circle-diag(diag(circle))
true_edge=edge_ind(circle,0.1^3)

rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

for(i in 1:100){
  samp=rmvnorm(n=100,sigma=circle)
  s=t(samp)%*%samp/100
  
  samp_edge=edge_ind(s,0.1^3)
  
  rmse_samp[i]=rmse(s,circle)
  mnorm_samp[i]=mnorm(s,circle)
  spec_samp[i]=spec_norm(s,circle)
  
  acc=accuracy(samp_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.2027399; 0.008246711

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.7727312; 0.1130067

#2norm
mean(spec_samp); sd(spec_samp)
#4.063132; 0.4248723

#specificity 
mean(sp_samp); sd(sp_samp)
#0.003659574; 0.001887596

#sensitivity
mean(se_samp); sd(se_samp)
# 1;0

#mcc
mean(mcc_samp); sd(mcc_samp)
#NaN; NA

#p=100, n=200
circle=2*diag(100)
for(i in 1:99){
  for(j in (i+1):100){
    if((j-i)==1){
      circle[i,j]=1
    }
    if(j==50){
      circle[1,j]=0.9
    }
  }
}
circle=t(circle)+circle-diag(diag(circle))
true_edge=edge_ind(circle,0.1^3)

rmse_samp=numeric(length=100); mnorm_samp=numeric(length=100); spec_samp=numeric(length=100)
sp_samp=numeric(length=100); se_samp=numeric(length=100); mcc_samp=numeric(length=100)

for(i in 1:100){
  samp=rmvnorm(n=200,sigma=circle)
  s=t(samp)%*%samp/200
  
  samp_edge=edge_ind(s,0.1^3)
  
  rmse_samp[i]=rmse(s,circle)
  mnorm_samp[i]=mnorm(s,circle)
  spec_samp[i]=spec_norm(s,circle)
  
  acc=accuracy(samp_edge,true_edge)
  sp_samp[i]=acc[1]
  se_samp[i]=acc[2]
  mcc_samp[i]=acc[3]
}

#rmse
mean(rmse_samp); sd(rmse_samp)
#0.1425057; 0.002964159

#mnorm
mean(mnorm_samp); sd(mnorm_samp)
#0.5878202; 0.07821426

#2norm
mean(spec_samp); sd(spec_samp)
# 4.267182; 0.2893713

#specificity 
mean(sp_samp); sd(sp_samp)
#0.005587629;  0.001111116

#sensitivity
mean(se_samp); sd(se_samp)
# 1;0

#mcc
mean(mcc_samp); sd(mcc_samp)
#NaN; NA
