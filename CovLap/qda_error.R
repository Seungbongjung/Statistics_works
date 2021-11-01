source("QDA.R")
source("model_choice&measure.R")
library(spcov)

parkinsons=as.matrix(read.csv("parkinsons.csv"))
true_status=as.numeric(parkinsons[,"status"])
parkinson=parkinsons[,-17]

#sample covariance
samp_error=qda_samp(true_status,parkinson,10)
mean(samp_error); sd(samp_error)
#0.1953846; 0.03424781

#graphical lasso
gl_error=qda_gl(true_status,parkinson,10)
mean(gl_error); sd(gl_error)
# 0.1961538; 0.0348284

#proposed 
proposed_error=qda_proposed(true_status,parkinson,10)
#

#MPP
mean(proposed_error[[1]]); sd(proposed_error[[1]])

#MAP
mean(proposed_error[[2]]); sd(proposed_error[[2]])