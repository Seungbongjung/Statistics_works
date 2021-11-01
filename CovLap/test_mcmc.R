source("mcmc_slabspike.R")

library(tictoc)
library(mvtnorm)
library(dplyr)

#true covariance
curr_ex=matrix(rep(0,12^2),nrow=12)
diag(curr_ex)=c(0.239,1.554,0.362,0.199,0.349,0.295,0.715,0.164,0.518,0.379,0.159,0.207)
curr_ex[1,2]=0.117; curr_ex[1,8]=0.031; curr_ex[3,4]=0.002; curr_ex[4,5]=0.094
curr_ex[5,12]=-0.036; curr_ex[6,7]=-0.229; curr_ex[6,8]=0.002; curr_ex[8,9]=0.112
curr_ex[8,10]=-0.028; curr_ex[8,11]=-0.008; curr_ex[9,10]=-0.193; curr_ex[9,11]=-0.090
curr_ex[10,11]=0.167

curr_ex=curr_ex+t(curr_ex)-diag(diag(curr_ex))

#sampling
s=rmvnorm(250,sigma=curr_ex)
s=t(s)%*%s/250

tic("sampling")
mcmc_samp=mcmc(250,s,1,250/2,0.1^3,10000,5000,curr_ex,2/11)
toc("done")

