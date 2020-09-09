##########################################################################
##### r code to show 2-stage regression and MLR numerically equivalent ####
###########################################################################
## step 1 simulate X from MVN 
set.seed(1234567)
n <- 100
R <- matrix(c(1,0.5,0.5,0.5,
              0.5,1,0.5,0.5,
              0.5,0.5,1,0.5,
              0.5,0.5,0.5,1), 
            nrow = 4, ncol = 4)

mu <- c(rep(0,4))
X<-mvtnorm::rmvnorm(n, mean = mu, sigma = R)
X ##### rectangle shape ###
### step 2: si,ulate Y from MNN
B<-c(runif(4,0,1))
dim(B)
Y<-mvtnorm::rmvnorm(n, mean = X%*%B, sigma = 0.03*diag(100))
dim(Y)
#### step 3: regress Y on X LS solution 
Beta<-solve((t(X)%*%X))%*%t(X)%*%Y
dim(Beta)
rowMeans(Beta)
#[1]  0.004407439 -0.006653397 -0.002431389  0.001406191
apply(Beta, 1, var)
#[1] 0.079245363 0.009194004 0.163940878 0.008255730


##### 2-stage linaer regression #####
beta_p<- function(p,X,Y) {
  Beta1<-solve((t(X[,-p])%*%X[,-p]))%*%t(X[,-p])%*%Y
  res1<-Y-X[,-p]%*%Beta1
  Beta2<-solve((t(X[,-p])%*%X[,-p]))%*%t(X[,-p])%*%X[,p]
  res2<-X[,p]-X[,-p]%*%Beta2
  Beta3<-solve((t(res2)%*%res2))%*%t(res2)%*%res1
  return(c(mean=mean(Beta3),
  var=apply(Beta3, 1, var)))
}
beta_p(1,X,Y)
#mean         var 
#0.004407439 0.079245363   
beta_p(2,X,Y)
#mean         var 
#-0.006653397  0.009194004  
beta_p(3,X,Y)
#mean        var 
#-0.002431389  0.163940878 
beta_p(4,X,Y)
#mean         var 
#0.001406191 0.008255730 
