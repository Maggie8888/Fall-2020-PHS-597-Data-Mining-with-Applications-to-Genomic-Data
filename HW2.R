##########################################################################################
###############  implement cyclic coordinate descent algorithm for lasso #################
##########################################################################################
### step 1 generate data ####
rm(list=ls())
n=100
set.seed(123456789)
X<-rnorm(n,2,1)
epsilon<-rnorm(n,0,0.01)
Y<-1+2*X+3*X^2+4*X^3+epsilon
predictors<-cbind(X,X^2,X^3,X^4,X^5,X^6,X^7,X^8,X^9,X^10)
xc = sweep(predictors, 2, colMeans(predictors))
sdc = sqrt(apply(xc, 2, crossprod)/nrow(predictors))
xs = sweep(xc, 2, sdc, "/")
ys = scale(Y)*sqrt(n/(n-1))

### step 2 define functions ####
cycliccoorddescent<-function(max.iter,x,y,beta_init,lambduh){
  beta<-matrix(rep(0,10*(max.iter+1)),nrow=max.iter+1,ncol=10)
  beta[1,]<-beta_init
  d = dim(x)[2]
  iter<- 2
  while (iter<=max.iter){
    for (j in 1:d){
      min_beta_j = min_beta_multivariate(x, y,beta[iter,],lambduh,j)
       beta[iter,j] = min_beta_j
    }
  #beta_vals = rbind(beta_init, beta)
  iter=iter+1
  beta[iter,] = beta[iter-1,]
  
  if (iter%%100 == 0) {
     print(paste0('Coordinate descent iteration', iter))
  }
  }  
  return (beta)
}

min_beta_multivariate<-function(x,y,betas,lambduh,j){
  n = length(y)
  norm_x_j = norm(as.matrix(x[,j]),"F")
  a = t(y- x[,-j]%*%betas[-j])%*%x[,j]
  passin = lambduh*n/2
  res<-soft_threshold(a,passin)
  return (res/(norm_x_j^2))
}

soft_threshold<-function(a,passin){
    if (-a > passin){
     return(a+passin)
  }
if (a > passin){
  return(a-passin)
}
    else {return(0)}
  
}

#### step 3 output ####

result<-cycliccoorddescent(1000,xs,ys,c(rep(0,10)),0.2)
result[1,]
result[1001,]

#[1] 0.0000000 0.1817000 0.7203633 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000


##########################################################################################
####################  compare results with glmnet ########################################
##########################################################################################
#install.packages("glmnet")
library(glmnet)
fit = glmnet(xs,ys,alpha = 1,intercept = F,standardize=F,lambda = rep(0.1,1000))
max(fit$lambda)
min(fit$lambda)
length(fit$lambda)
fit$beta[,1000]

#0.0000000 0.1817000 0.7203633 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 

