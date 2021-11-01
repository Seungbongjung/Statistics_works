#protein data
library(sparsebn)

dat=sparsebnData(cytometryContinuous$data, type = "c", ivn =NULL)
true_cov=estimate.covariance(dat)
