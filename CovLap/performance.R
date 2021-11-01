# edge count
edge_ind=function(m,tol){
  edge_list=m[upper.tri(m)]
  edge_indicator=as.numeric((abs(edge_list)>=tol))
  return(edge_indicator)
}

# performance check
performance1=function(esti,m,tol){
  edge_list_est=esti[upper.tri(esti)]
  edge_list_m=m[upper.tri(m)]
  
  edge_ind=as.numeric(abs(edge_list_est)>=tol)
  true_edge=as.numeric(abs(edge_list_m)>=tol)
  
  edge_ind1=which(edge_ind==1)
  edge_ind0=which(edge_ind==0)
  
  true_edge1=which(true_edge==1)
  true_edge0=which(true_edge==0)
  
  
  tp=length(intersect(edge_ind1,true_edge1))
  tn=length(intersect(edge_ind0,true_edge0))
  fp=length(intersect(edge_ind1,true_edge0))
  fn=length(intersect(edge_ind0,true_edge1))
  
  sp=tn/(tn+fp)
  se=tp/(tp+fn)
  return(c(sp,se))
}

performance2=function(esti,cri,tol){
  edge_list_est=esti[upper.tri(esti)]
  
  edge_ind=as.numeric(abs(edge_list_est)>=tol)
  
  edge_ind1=which(edge_ind==1)
  edge_ind0=which(edge_ind==0)
  
  true_edge1=which(cri==1)
  true_edge0=which(cri==0)
  
  
  tp=length(intersect(edge_ind1,true_edge1))
  tn=length(intersect(edge_ind0,true_edge0))
  fp=length(intersect(edge_ind1,true_edge0))
  fn=length(intersect(edge_ind0,true_edge1))
  
  sp=tn/(tn+fp)
  se=tp/(tp+fn)
  return(c(sp,se))
}