source("bcd_algorithm.R")
source("performance.R")
library(mvtnorm)

#sparse covariance that mimics daily currency exchange rate structure 
#p=12, number of nodes=13

curr_ex=matrix(rep(0,12^2),nrow=12)
diag(curr_ex)=c(0.239,1.554,0.362,0.199,0.349,0.295,0.715,0.164,0.518,0.379,0.159,0.207)
curr_ex[1,2]=0.117; curr_ex[1,8]=0.031; curr_ex[3,4]=0.002; curr_ex[4,5]=0.094
curr_ex[5,12]=-0.036; curr_ex[6,7]=-0.229; curr_ex[6,8]=0.002; curr_ex[8,9]=0.112
curr_ex[8,10]=-0.028; curr_ex[8,11]=-0.008; curr_ex[9,10]=-0.193; curr_ex[9,11]=-0.090
curr_ex[10,11]=0.167

curr_ex=curr_ex+t(curr_ex)-diag(diag(curr_ex))

#sample from true covariance, n=250
set.seed(20210706)
samp=rmvnorm(250,sigma=curr_ex)
s=t(samp)%*%samp/250

#edge of true covariance
edge0=edge_ind(curr_ex,0.1^3)

r1=sample(1:choose(12,2),replace=FALSE,size=2)
edge1=edge0; edge1[r1]=1-edge1[r1]

r2=sample(1:choose(12,2),replace=FALSE,size=3)
edge2=edge0; edge2[r2]=1-edge2[r2]

r3=sample(1:choose(12,2),replace=FALSE,size=4)
edge3=edge0; edge3[r3]=1-edge1[r3]

r4=sample(1:choose(12,2),replace=FALSE,size=5)
edge4=edge0; edge4[r4]=1-edge4[r4]

#performance

#edge0

#bcd1
cov1_edge0=bcd_covariance1(s,250,0.005,1,250/2,10000,0.1^6,edge0)

#bcd2
cov2_edge0=bcd_covariance2(s,250,1,250/2,10000,0.1^6,edge0)

#bcd1 performance
cov1_edge0_performance=performance(cov1_edge0,curr_ex,0.1^3)
cov1_edge0_performance
#0.4528302 0.9230769

#bcd2 performance
cov2_edge0_performance=performance(cov2_edge0,curr_ex,0.1^3)
cov2_edge0_performance
#1.0000000 0.9230769

#edge1

#bcd1
cov1_edge1=bcd_covariance1(s,250,0.005,1,250/2,10000,0.1^3,edge1)

#bcd2
cov2_edge1=bcd_covariance2(s,250,1,250/2,10000,0.1^3,edge1)

#bcd1 performance
cov1_edge1_performance=performance2(cov1_edge1,edge1,0.1^3)
cov1_edge1_performance
#0.4313725 0.9333333

#bcd2 performance
cov2_edge1_performance=performance2(cov2_edge1,edge1,0.1^3)
cov2_edge1_performance
#1 1

#edge2

#bcd1
cov1_edge2=bcd_covariance1(s,250,0.005,1,250/2,10000,0.1^3,edge2)

#bcd2
cov2_edge2=bcd_covariance2(s,250,1,250/2,10000,0.1^3,edge2)

#bcd1 performance
cov1_edge2_performance=performance2(cov1_edge2,edge2,0.1^3)
cov1_edge2_performance
#0.4423077 0.9285714

#bcd2 performance
cov2_edge2_performance=performance2(cov2_edge2,edge2,0.1^3)
cov2_edge2_performance
#1.0000000 0.8571429

#edge3

#bcd1
cov1_edge3=bcd_covariance1(s,250,0.005,1,250/2,10000,0.1^3,edge3)

#bcd2
cov2_edge3=bcd_covariance2(s,250,1,250/2,10000,0.1^3,edge3)

#bcd1 performance
cov1_edge3_performance=performance2(cov1_edge3,edge3,0.1^3)
cov1_edge3_performance
#0.4528302 0.9230769

#bcd2 performance
cov2_edge3_performance=performance2(cov2_edge3,edge3,0.1^3)
cov2_edge3_performance
#1.0000000 0.9230769
