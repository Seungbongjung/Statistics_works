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

#proposed
sssl_perf1=bayesian_measure(200,1,200*0.06,0.1^3,12000,3000,4/99,1,toep)
sssl_perf1[[1]];sssl_perf1[[2]];sssl_perf1[[3]];sssl_perf1[[4]];sssl_perf1[[5]];sssl_perf1[[6]]
sssl_perf1[[7]];sssl_perf1[[8]];sssl_perf1[[9]];sssl_perf1[[10]];sssl_perf1[[11]];sssl_perf1[[12]]

#MPP
#rmse
mean(sssl_perf1[[1]]); sd(sssl_perf1[[1]])
#0.08556214

#mnorm
mean(sssl_perf1[[2]]); sd(sssl_perf1[[2]])
#0.4973146

#2norm
mean(sssl_perf1[[3]]); sd(sssl_perf1[[3]])
#1.683308

#specificity
mean(sssl_perf1[[4]]); sd(sssl_perf1[[4]])
#  0.9979403

#sensitivity
mean(sssl_perf1[[5]]); sd(sssl_perf1[[5]])
#0.1259843

#mcc
mean(sssl_perf1[[6]]); sd(sssl_perf1[[6]])
#NaN; NA

#MAP
#rmse
mean(sssl_perf1[[7]]); sd(sssl_perf1[[7]])
#0.08448291

#mnorm
mean(sssl_perf1[[8]]); sd(sssl_perf1[[8]])
#0.5308378

#2norm
mean(sssl_perf1[[9]]); sd(sssl_perf1[[9]])
#1.699444

#specificity
mean(sssl_perf1[[10]]); sd(sssl_perf1[[10]])
#0.9752832

#sensitivity
mean(sssl_perf1[[11]]); sd(sssl_perf1[[11]])
#0.1732283

#mcc
mean(sssl_perf1[[12]]); sd(sssl_perf1[[12]])
#NaN; NA

#scenario6
model=read.csv("model.csv")
model=model[,-1]
model=as.matrix(model)

sssl_perf2=bayesian_measure(100,1,100*0.06,0.1^3,12000,3000,4/49,1,model)
sssl_perf2[[1]];sssl_perf2[[2]];sssl_perf2[[3]];sssl_perf2[[4]];sssl_perf2[[5]];sssl_perf2[[6]]
sssl_perf2[[7]];sssl_perf2[[8]];sssl_perf2[[9]];sssl_perf2[[10]];sssl_perf2[[11]];sssl_perf2[[12]]

#MPP
#rmse
mean(sssl_perf2[[1]]); sd(sssl_perf2[[1]])
#0.0814121

#mnorm
mean(sssl_perf2[[2]]); sd(sssl_perf2[[2]])
#0.3222389

#2norm
mean(sssl_perf2[[3]]); sd(sssl_perf2[[3]])
#1.873908

#specificity
mean(sssl_perf2[[4]]); sd(sssl_perf2[[4]])
#0.9979403

#sensitivity
mean(sssl_perf2[[5]]); sd(sssl_perf2[[5]])
#0.08267717

#mcc
mean(sssl_perf2[[6]]); sd(sssl_perf2[[6]])
#NaN; NA

#MAP
#rmse
mean(sssl_perf2[[7]]); sd(sssl_perf2[[7]])
#0.08200986

#mnorm
mean(sssl_perf2[[8]]); sd(sssl_perf2[[8]])
#0.3944104

#2norm
mean(sssl_perf2[[9]]); sd(sssl_perf2[[9]])
#1.788268

#specificity
mean(sssl_perf2[[10]]); sd(sssl_perf2[[10]])
#0.9783728

#sensitivity
mean(sssl_perf2[[11]]); sd(sssl_perf2[[11]])
#0.1259843

#mcc
mean(sssl_perf2[[12]]); sd(sssl_perf2[[12]])
#NaN; NA