#scenario3

#p=50
random1=matrix(rep(0,50^2),nrow=50)
random1_off=rep(0,choose(50,2))

for(i in 1:choose(50,2)){
  samp=sample(c(0,1),size=1,replace=FALSE,prob=c(0.98,0.02))
  if(samp==1){
    random1_off[i]=sample(c(-1,1),size=1,replace=FALSE,prob=c(0.5,0.5))
  }
}

random1[upper.tri(random1)]=random1_off
random1=random1+t(random1)
x=eigen(random1)$values
random1=random1+diag(rep(2.536,50))
write.csv(random1,"random_structure2_50.csv")

#p=100
random1=matrix(rep(0,100^2),nrow=100)
random1_off=rep(0,choose(100,2))

for(i in 1:choose(100,2)){
  samp=sample(c(0,1),size=1,replace=FALSE,prob=c(0.98,0.02))
  if(samp==1){
    random1_off[i]=sample(c(-1,1),size=1,replace=FALSE,prob=c(0.5,0.5))
  }
}

random1[upper.tri(random1)]=random1_off
random1=random1+t(random1)

random1=random1+diag(rep(3.111,100))
write.csv(random1,"random_structure2_100.csv")

#scenario4
#p=50
random2=matrix(rep(0,50^2),nrow=50)

for(i in 1:49){
  random2[i,i+1]=0.4
}
random2=random2+t(random2)  

random2=random2+diag(rep(0.831,50))
write.csv(random2,"random_structure3_50.csv")

#p=100
random2=matrix(rep(0,100^2),nrow=100)
               
for(i in 1:99){
  random2[i,i+1]=0.4
}
random2=random2+t(random2)  
random2=random2+diag(rep(0.8157,100))
write.csv(random2,"random_structure3_100.csv")
eigen.spam(random2,nev=50)

#scenario2 

#p=50
random3=matrix(rep(0,50^2),nrow=50)
random3_off=rep(0,choose(50,2))
sparse=sample(1:choose(50,2),size=choose(50,2)*0.2,replace=FALSE)
off_diag=runif(choose(50,2)*0.2,min=0,max=0.02)
random3_off[sparse]=off_diag
random3[upper.tri(random3)]=random3_off
random3=random3+t(random3)
random3=random3+diag(rgamma(50,1,1))
write.csv(random3,"random_structure1_50.csv")

#p=100
random3=matrix(rep(0,100^2),nrow=100)
random3_off=rep(0,choose(100,2))
sparse=sample(1:choose(100,2),size=choose(100,2)*0.2,replace=FALSE)
off_diag=runif(choose(100,2)*0.2,min=0,max=0.02)
random3_off[sparse]=off_diag
random3[upper.tri(random3)]=random3_off
random3=random3+t(random3)
random3=random3+diag(rgamma(100,1,1))
write.csv(random3,"random_structure1_100.csv")
