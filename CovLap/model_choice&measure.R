library(dplyr)
library(mvtnorm)
source("mcmc_slabspike.R")

#median probability model
mpp=function(mat){
  r=nrow(mat); c=ncol(mat)
  result=numeric(length=c)
  for(i in 1:c){
    result[i]=as.numeric(mean(mat[,i])>=0.5)
  }
  return(result)
}

#maximum a probability
mmp=function(mat){

  temp=as.data.frame(mat)
  temp=temp%>%group_by_all%>%count
  l=ncol(temp)
  
  n=temp$n
  ind=which(n==max(n))
  model=as.numeric(temp[ind,1:(l-1)])
  return(model)
}

#rmse
rmse=function(mat1,mat2){
  p=nrow(mat1)
  val=sqrt(sum((mat1-mat2)^2))/p
  return(val)
}

#mnorm
mnorm=function(mat1,mat2){
  temp=abs(mat1-mat2)
  val=max(temp)
  return(val)
}

#spectral norm 
spec_norm=function(mat1,mat2){
  temp=mat1-mat2
  temp=t(temp)%*%temp
  val=sqrt(max(eigen(temp)$values))
  return(val)
}

#accuracy
accuracy=function(edge1,edge2){
  
  edge_ind1=which(edge1==1)
  edge_ind0=which(edge1==0)
  
  true_edge1=which(edge2==1)
  true_edge0=which(edge2==0)
  
  
  tp=length(intersect(edge_ind1,true_edge1))
  tn=length(intersect(edge_ind0,true_edge0))
  fp=length(intersect(edge_ind1,true_edge0))
  fn=length(intersect(edge_ind0,true_edge1))
  
  sp=tn/(tn+fp)
  se=tp/(tp+fn)
  mcc=(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  return(c(sp,se,mcc))
  
}

#performance 
bayesian_measure=function(n,v,lambda,tol,size,burnin,q,rep,true_cov){
  
  temp1=numeric(length=rep); temp2=numeric(length=rep)
  temp3=numeric(length=rep); temp4=numeric(length=rep)
  temp5=numeric(length=rep); temp6=numeric(length=rep)
  temp7=numeric(length=rep); temp8=numeric(length=rep)
  temp9=numeric(length=rep); temp10=numeric(length=rep)
  temp11=numeric(length=rep); temp12=numeric(length=rep)

  true_edge=edge_ind(true_cov,0.1^3)

  for(i in 1:rep){
   print(paste(c("replication: ",i),sep=""))
   samp=rmvnorm(n,sigma=true_cov)
   s=t(samp)%*%samp/n
   mcmc_samp=mcmc(n,s,v,lambda,tol,size,burnin,q)
   mcmc_mpp=mpp(mcmc_samp)
   mcmc_mmp=mmp(mcmc_samp)
   mpp_bcd=bcd_covariance2(s,n,v,lambda,10000,tol,mcmc_mpp)
   mmp_bcd=bcd_covariance2(s,n,v,lambda,10000,tol,mcmc_mmp)
   
   temp1[i]=rmse(mpp_bcd,true_cov)
   temp2[i]=mnorm(mpp_bcd,true_cov)
   temp3[i]=spec_norm(mpp_bcd,true_cov)
   
   acc_mpp=accuracy(mcmc_mpp,true_edge)
   
   temp4[i]=acc_mpp[1]
   temp5[i]=acc_mpp[2]
   temp6[i]=acc_mpp[3]
   
   temp7[i]=rmse(mmp_bcd,true_cov)
   temp8[i]=mnorm(mmp_bcd,true_cov)
   temp9[i]=spec_norm(mmp_bcd,true_cov)
   
   acc_mmp=accuracy(mcmc_mmp,true_edge)
   
   temp10[i]=acc_mmp[1]
   temp11[i]=acc_mmp[2]
   temp12[i]=acc_mmp[3]
  }
  measure=list(mpp_rmse=temp1,mpp_mnorm=temp2,mpp_spec=temp3,mpp_sp=temp4,mpp_se=temp5,
               mpp_mcc=temp6, mmp_rmse=temp7, mmp_mnorm=temp8, mmp_spec=temp9, mmp_sp=temp10,
               mmp_se=temp11,mmp_mcc=temp12)
  return(measure)
}

