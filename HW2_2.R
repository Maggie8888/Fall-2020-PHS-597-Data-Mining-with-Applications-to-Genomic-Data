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
xc = sweep(predictors, 2, colMeans(predictors))
sdc = sqrt(apply(xc, 2, crossprod)/nrow(predictors))
xs = sweep(xc, 2, sdc, "/")
ys = scale(Y)*sqrt(n/(n-1))

### step 2 define functions ##########
##PLS   Partial Least Squares Regrassion
##between the independent variables, X and dependent Y as
## X = T*P' + E;
## Y = U*Q' + F = T*B*Q' + F1;
##
## Inputs:
## X     data matrix of independent variables
## Y     data matrix of dependent variables
## tol1  the tolerant of convergence 
## tol2  the tolerant of convergence 
###
## Outputs:
## T     score matrix of X
## P     loading matrix of X
## U     score matrix of Y
## Q     loading matrix of Y
## B     matrix of regression coefficient
## W     weight matrix of X
## Using the PLS model, for new X1, Y1 can be predicted as
## Y1 = (X1*P)*B*Q' = X1*(P*B*Q')
## or
## Y1 = X1*(W*inv(P'*W)*inv(T'*T)*T'*Y)
##
## Without Y provided, the function will return the principal components as
## X = T*P' + E
##
PLS<-function(X,Y,tol1,tol2){
  rX<-dim(X)[1]
  cX<-dim(X)[2]
  rY<-dim(Y)[1]
  cY<-dim(Y)[2]
## Allocate memory to the maximum size 
n=max(cX,cY)
T=matrix(rep(0,rX*n),nrow=rX,ncol=n)
P=matrix(rep(0,cX*n),nrow=cX,ncol=n)
U=matrix(rep(0,rY*n),nrow=rY,ncol=n)
Q=matrix(rep(0,r=cY*n),nrow=cY,ncol=n)
B=matrix(rep(0,n*n),nrow=n,ncol=n)
W=P
k=0
## iteration loop if residual is larger than specfied
while (norm(Y)>tol1 && k<n && !is.na(norm(Y,"2"))){
## choose the column of x has the largest square of sum as t.
## choose the column of y has the largest square of sum as u.    
tidx =  which.max(apply(X, 2, function(x) sum(t(x)%*%x)))
uidx =   which.max(apply(Y, 2, function(y) sum(t(y)%*%y)))
t1 = X[,tidx]
u = Y[,uidx]
t = matrix(rep(0,rX*1),nrow=rX,ncol=1)

## iteration for outer modeling until convergence
while (norm(t1-t) > tol2){
  w = t(X)%*%u
  w = w/norm(w,"2")
  t = t1
  t1 = X%*%w;
  q = t(Y)%*%t1
q = q/norm(q);
u = Y%*%q
}


##update p based on t
t=t1
p=t(X)%*%t
p=p/(rep(t(t)%*%t,length(p)))
pnorm=norm(p,"2")
p=p/pnorm
t=pnorm*t
w=pnorm*w

## regression and residuals
b = (t(u)%*%t)/(t(t)%*%t)
X = X - t%*%t(p)
Y = Y - as.numeric(b)*t%*%t(q)
## save iteration results to outputs:
  k=k+1
  T[,k]=t
  P[,k]=p
  U[,k]=u
  Q[,k]=q
  W[,k]=w
  B[k,k]=b

}
  #B=B[1:k,1:k]
  return(list(T=T,P=P,U=U,Q=Q,B=B,W=W))
}

res<-PLS(xs,ys,0.0000001,0.0000001)
res$B
res$W%*%solve(t(res$P)%*%res$W)%*%t(res$Q)
#install.packages("pls")
library(pls)
mod <- plsr(ys ~ xs)
mod$coefficients
