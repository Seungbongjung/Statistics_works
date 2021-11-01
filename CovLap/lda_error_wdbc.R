source("lda_cv_wdbc.R")
source("model_choice&measure.R")
library(spcov)

#dataset

wdbc=read.csv("wdbc.csv",header=FALSE)
true_status=wdbc[,1]
true_status[which(true_status=="M")]=1; true_status[which(true_status=="B")]=0
true_status=as.numeric(true_status)
wdbc=as.matrix(wdbc[,-1])

#error 
error=lda_cv(true_status,wdbc,10)

#sample covariance matrix
mean(error[[1]]); sd(error[[1]])
#0.07698413; 0.01933812

#graphical lasso
mean(error[[2]]); sd(error[[2]])
# 0.07222222; 0.01670605

#proposed (MPP)
mean(error[[3]]); sd(error[[3]])
#0.06613757; 0.01635521

#proposed (MAP)
mean(error[[4]]); sd(error[[4]])
#0.06613757; 0.01635521
