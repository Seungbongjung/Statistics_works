source("model_choice&measure_ver2.R")

#scenario1
#sparse covariance that mimics currency exchange rate, n=250
curr_ex=matrix(rep(0,12^2),nrow=12)
diag(curr_ex)=c(0.239,1.554,0.362,0.199,0.349,0.295,0.715,0.164,0.518,0.379,0.159,0.207)
curr_ex[1,2]=0.117; curr_ex[1,8]=0.031; curr_ex[3,4]=0.002; curr_ex[4,5]=0.094
curr_ex[5,12]=-0.036; curr_ex[6,7]=-0.229; curr_ex[6,8]=0.002; curr_ex[8,9]=0.112
curr_ex[8,10]=-0.028; curr_ex[8,11]=-0.008; curr_ex[9,10]=-0.193; curr_ex[9,11]=-0.090
curr_ex[10,11]=0.167
curr_ex=t(curr_ex)+curr_ex-diag(diag(curr_ex))

set.seed(1)
sssl_perf=bayesian_measure(250,1,250*0.06,0.1^3,12000,3000,4/11,100,curr_ex)

#median probability model
#rmse
mean(sssl_perf[[1]]); sd(sssl_perf[[1]])
#0.01992495; 0.006090275
#0.02032584l; 0.005963799

#mnorm
mean(sssl_perf[[2]]); sd(sssl_perf[[2]])
#0.1473022; 0.07214723
#0.1535876; 0.08045469

#2norm
mean(sssl_perf[[3]]);sd(sssl_perf[[3]])
#0.1781454; 0.07596491
#0.1806409; 0.079999

#specificity
mean(sssl_perf[[4]]); sd(sssl_perf[[4]])
#0.9986792; 0.004838347
#0.9981132; 0.005688893

#sensitivity
mean(sssl_perf[[5]]); sd(sssl_perf[[5]])
#0.5515385; 0.06473349
#0.5569231; 0.06386436
#mcc
mean(sssl_perf[[6]]); sd(sssl_perf[[6]])
#0.6997885; 0.0478887
#0.7018462; 0.04833532

#maximum a posteriori (MAP)
#rmse
mean(sssl_perf[[7]]); sd(sssl_perf[[7]])
#0.01995951;  0.006098553
#0.0201802l;  0.00608459
#mnorm
mean(sssl_perf[[8]]); sd(sssl_perf[[8]])
#0.1468436; 0.07248782
#0.153361; 0.08114199

#2norm
mean(sssl_perf[[9]]);sd(sssl_perf[[9]])
#0.1779943; 0.07602901
#0.1787265; 0.08045904

#specificity
mean(sssl_perf[[10]]); sd(sssl_perf[[10]])
#0.9984906; 0.005144527
#0.9971698; 0.006771137

#sensitivity
mean(sssl_perf[[11]]); sd(sssl_perf[[11]])
#0.5515385; 0.06473349
#0.5638462; 0.06744658

#mcc
mean(sssl_perf[[12]]); sd(sssl_perf[[12]])
#0.6991698; 0.04829839
#0.7036343; 0.05295202

#scenario2
#random structure1-dimension p=50, n=100, mu=0.05
random1=read.csv("random_structure1.csv")
random1=as.matrix(random1[,-1])

sssl_perf1=bayesian_measure(100,1,100*0.06,0.1^3,12000,3000,4/49,1,random1)
sssl_perf1[[1]];sssl_perf1[[2]];sssl_perf1[[3]];sssl_perf1[[4]];sssl_perf1[[5]];sssl_perf1[[6]]
sssl_perf1[[7]];sssl_perf1[[8]];sssl_perf1[[9]];sssl_perf1[[10]];sssl_perf1[[11]];sssl_perf1[[12]]


#median probability model
#rmse
mean(sssl_perf1[[1]]); sd(sssl_perf1[[1]])
#0.05120519
#0.04587659

#mnorm
mean(sssl_perf1[[2]]); sd(sssl_perf1[[2]])
#0.9823065
#1.159046

#2norm
mean(sssl_perf1[[3]]);sd(sssl_perf1[[3]])
#1.477123
#1.190308

#specificity
mean(sssl_perf1[[4]]); sd(sssl_perf1[[4]])
#0.9116751
#0.9949239

#sensitivity
mean(sssl_perf1[[5]]); sd(sssl_perf1[[5]])
#0.1083333
#0.03333333

#mcc
mean(sssl_perf1[[6]]); sd(sssl_perf1[[6]])
#NaN; NA
#NaN; NA

#maximum a posteriori (MAP)
#rmse
mean(sssl_perf1[[7]]); sd(sssl_perf1[[7]])
#0.06603483
#0.04634315

#mnorm
mean(sssl_perf1[[8]]); sd(sssl_perf1[[8]])
#1.011064
#1.158608

#2norm
mean(sssl_perf1[[9]]);sd(sssl_perf1[[9]])
#1.415001
#1.185584

#specificity
mean(sssl_perf1[[10]]); sd(sssl_perf1[[10]])
#0.742132
#0.986802

#sensitivity
mean(sssl_perf1[[11]]); sd(sssl_perf1[[11]])
#0.2333333
#0.0375

#mcc
mean(sssl_perf1[[12]]); sd(sssl_perf1[[12]])
#NaN; NA

#scenario3
#random structure2-dimension p=100, n=200
random2=read.csv("random_structure2.csv",header=FALSE)
random2=as.matrix(random2)

set.seed(3)
library(tictoc)
tic()
sssl_perf2=bayesian_measure(200,1,200*0.06,1,12000,3000,1/2,1,random2)
toc()
sssl_perf2[[1]];sssl_perf2[[2]];sssl_perf2[[3]];sssl_perf2[[4]];sssl_perf2[[5]];sssl_perf2[[6]]
sssl_perf2[[7]];sssl_perf2[[8]];sssl_perf2[[9]];sssl_perf2[[10]];sssl_perf2[[11]];sssl_perf2[[12]]

#median probability model
#rmse
mean(sssl_perf1[[1]]); sd(sssl_perf1[[1]])
#0.05120519

#mnorm
mean(sssl_perf1[[2]]); sd(sssl_perf1[[2]])
#0.9823065

#2norm
mean(sssl_perf1[[3]]);sd(sssl_perf1[[3]])
#1.477123

#specificity
mean(sssl_perf1[[4]]); sd(sssl_perf1[[4]])
#0.9116751

#sensitivity
mean(sssl_perf1[[5]]); sd(sssl_perf1[[5]])
#0.1083333

#mcc
mean(sssl_perf1[[6]]); sd(sssl_perf1[[6]])
#NaN; NA

#maximum a posteriori (MAP)
#rmse
mean(sssl_perf1[[7]]); sd(sssl_perf1[[7]])
#0.06603483

#mnorm
mean(sssl_perf1[[8]]); sd(sssl_perf1[[8]])
#1.011064

#2norm
mean(sssl_perf1[[9]]);sd(sssl_perf1[[9]])
#1.415001

#specificity
mean(sssl_perf1[[10]]); sd(sssl_perf1[[10]])
#0.742132

#sensitivity
mean(sssl_perf1[[11]]); sd(sssl_perf1[[11]])
#0.2333333

#mcc
mean(sssl_perf1[[12]]); sd(sssl_perf1[[12]])
#NaN; NA

#scenario4
#random structure3-dimension p=100, n=200
random3=read.csv("random_structure3.csv",header=FALSE)
random3=as.matrix(random3)
set.seed(4)
tic()
sssl_perf3=bayesian_measure(200,1,200*0.06,1,1,0,4/99,1,random3)
toc()
sssl_perf3[[1]];sssl_perf3[[2]];sssl_perf3[[3]];sssl_perf3[[4]];sssl_perf3[[5]];sssl_perf3[[6]]
sssl_perf3[[7]];sssl_perf3[[8]];sssl_perf3[[9]];sssl_perf3[[10]];sssl_perf3[[11]];sssl_perf3[[12]]


#real dataset1: protein data

#real dataset2: 

