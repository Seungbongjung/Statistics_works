source("LDA.R")
source("model_choice&measure.R")

parkinsons=as.matrix(read.csv("parkinsons.csv"))
true_status=as.numeric(parkinsons[,"status"])
parkinson=parkinsons[,-17]

#sample covariance matrix
samp_error=lda_samp(true_status,parkinson)
samp_error
#0.1282051

#graphical lasso
gl_error=lda_gl(true_status,parkinson)
gl_error
#0.1384615

#proposed
proposed_error=lda_proposed(true_status,parkinson)


#MPP
proposed_error[1]

#MAP
proposed_error[2]

