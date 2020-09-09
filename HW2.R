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
y<-Y-mean(Y)
x<-matrix(rep(0,1000),ncol=10,nrow=100)
for (j in 1:10){
  x[1:100,j]<-(predictors[1:100,j]-mean(predictors[1:100,j]))/sd(predictors[1:100,j])
}

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

result<-cycliccoorddescent(1000,x,y,c(rep(0,10)),0.2)
result[1,]
result[1001,]

#[1]  0.00000 22.26865 37.91739  8.15412  0.00000  0.00000  0.00000  0.00000  0.00000  0.00000

##########################################################################################
####################  compare results with glmnet ########################################
##########################################################################################
#install.packages("glmnet")
library(glmnet)
fit1 = glmnet(x,y,alpha = 1,intercept = F,standardize=F,lambda = rep(0.2,270))
max(fit1$lambda)
min(fit1$lambda)
length(fit1$lambda)
fit1$beta[,270]

#V1        V2        V3        V4        V5        V6        V7        V8        V9       V10 
#0.000000 22.208013 37.953496  8.074919  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000 

