##########################################################################
##### r code to show 2-stage regression and MLR numerically equivalent ####
###########################################################################
## step 1 simulate X from MVN 
rm(list=ls())
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
### step 2: simulate Y from MNN
B<-c(runif(4,0,1))
dim(B)
sigma<-rnorm(100)
Y<-X%*%B+sigma
dim(Y)
Mydata<-as.data.frame(cbind(Y,X))
lm <- lm(Y~X[,1]+X[,2]+X[,3]+X[,4],data=Mydata) 
sum<-summary(lm)
uci.lm<-sum$coefficients[,1]+2*sum$coefficients[,2]
lci.lm<-sum$coefficients[,1]-2*sum$coefficients[,2]
res.lm<-rbind(lci.lm,uci.lm)
res.lm
#(Intercept)      X[, 1]    X[, 2]     X[, 3]    X[, 4]
#lci.lm  -0.1283217 0.009503603 0.2506026 -0.2748928 0.6265862
#uci.lm   0.2556606 0.496723165 0.7345665  0.2889247 1.0687281
### CI of the coefficients from regression ####

#### parametric bootstrap CI 

#sampling rows from the data frame:
resample <- function(x) {sample(x,size=length(x),replace=TRUE)}
resample.Mydata <- function(Mydata) {
  sample.rows <- resample(1:nrow(Mydata))
  return(sample.rows)
}


est<- function(subset,data=Mydata) 
{ fit <- lm(Y ~ X[,1]+X[,2]+X[,3]+X[,4],data=data,subset=subset) 
return(coefficients(fit))
}

cis.par <- function(B,alpha) 
  { tboot <- replicate(B,est(resample.Mydata(Mydata))) 
  low.quantiles <- apply(tboot,1,quantile,probs=alpha/2) 
  high.quantiles <- apply(tboot,1,quantile,probs=1-alpha/2) 
  low.cis <- 2*coefficients(lm) - high.quantiles 
  high.cis <- 2*coefficients(lm) - low.quantiles
cis <- rbind(low.cis,high.cis)
rownames(cis) <- as.character(c(alpha/2,1-alpha/2)) 
return(cis)
}

set.seed(1234)
signif(cis.par(B=1e4,alpha=0.05),3)

#      (Intercept) X[, 1] X[, 2] X[, 3] X[, 4]
#0.025      -0.129 0.0244  0.279 -0.267  0.648
#0.975       0.257 0.5040  0.720  0.280  1.050

###### non-parametric bootstrap 
### using LS estimation ####

beta_p<- function(p,X,Y) {
  Beta1<-ginv((t(X[,-p])%*%X[,-p]))%*%t(X[,-p])%*%Y
  res1<-Y-X[,-p]%*%Beta1
  Beta2<-ginv((t(X[,-p])%*%X[,-p]))%*%t(X[,-p])%*%X[,p]
  res2<-X[,p]-X[,-p]%*%Beta2
  Beta3<-ginv((t(res2)%*%res2))%*%t(res2)%*%res1
  return(mean=mean(Beta3))
}

rege.data<- function(Mydata){ 
  subset<-Mydata[resample.Mydata(Mydata),]
  Y<-subset[,1]
  X<-subset[,2:5]
  X<-as.matrix(X)
  return(list(X=X,Y=Y))
}

npr.cis<-function(B,alpha,res){ 
  tboot<-matrix(NA,nrow=B,ncol=4)
  tboot[,1] <- replicate(B,beta_p(1,rege.data(Mydata)$X,rege.data(Mydata)$Y)) 
  tboot[,2] <- replicate(B,beta_p(2,rege.data(Mydata)$X,rege.data(Mydata)$Y))
  tboot[,3] <- replicate(B,beta_p(3,rege.data(Mydata)$X,rege.data(Mydata)$Y))
  tboot[,4] <- replicate(B,beta_p(4,rege.data(Mydata)$X,rege.data(Mydata)$Y))
  low.quantiles <- apply(tboot,2,quantile,probs=alpha/2) 
  high.quantiles <- apply(tboot,2,quantile,probs=1-alpha/2) 
  low.cis <- 2*res - high.quantiles 
  high.cis <- 2*res - low.quantiles
  cis <- rbind(low.cis,high.cis) 
  return(list(cis=cis,tboot=t(tboot)))
}
res<-beta_p(1,X,Y)
set.seed(1234)
res.cis <-npr.cis(1e4,0.05,res)
res.cis$cis
#### wider CI #####
#                 [,1]      [,2]        [,3]      [,4]
#low.cis  0.0114205 0.0287247 -0.02802851 0.0723781
#high.cis 0.9200968 0.9218671  0.99537253 0.8871391
