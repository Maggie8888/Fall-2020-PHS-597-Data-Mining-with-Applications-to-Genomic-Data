##########################################################################
################         spline model in R            ####################
###########################################################################
## step 1 simulate data set 
rm(list=ls())
set.seed(123456)
x<-rnorm(100,0,1)
epsilon<-rnorm(100,0,1)
y<-exp(x)/(1+exp(x))+epsilon

# define new predictors
x0 <- (x - 0)
x0[x0<0]<- 0
x1 <-(x+1)
x1[x1<0] <- 0
x2 <- (x - 1)
x2[x2<0]<- 0

#1) cubic splines;knots at 0, -1 and 1
x0.cubed <- x0^3
x1.cubed <- x1^3
x2.cubed <- x2^3
x.squared <- x^2
x.cubed <- x^3
fit <- lm( y ~ x + x.squared + x.cubed + x0.cubed + x1.cubed + x2.cubed )
print(summary(fit))
#### LS methods for coefficients ####
X<-cbind(rep(1,100),x,x.squared,x.cubed,x0.cubed, x1.cubed, x2.cubed)

Beta<-solve((t(X)%*%X))%*%t(X)%*%y
dim(Beta)
rowMeans(Beta)
#x   x.squared     x.cubed    x0.cubed    x1.cubed    x2.cubed 
#1.47516834  3.77236823  3.54377091  0.84611200  0.04849454 -1.14287354  1.02233256 
#### comparing with package 
#install.packages("splines")


#2) natural cubic splines;knots at 0, -1 and 1
X<-cbind(rep(1,100), x,((x1.cubed-x2.cubed)/2-(x0.cubed-x2.cubed)/1))
Beta<-solve((t(X)%*%X))%*%t(X)%*%y
dim(Beta)
rowMeans(Beta)
## fit linear model in expanded variables
x2<-((x0.cubed-x2.cubed)/1-(x1.cubed-x2.cubed)/2)
fit_e  <- lm(y ~x+x2)
## get estimate of beta
print(beta_e <- coef(fit_e))
#### comparing with package 
summary(lm(y ~ ns(x, df = 3,knots = c(-1,0,1))))
#3) smoothing splines;knots at 0, -1 and 1

smooth.matrix = function(x, df){
  ## return the smoother matrix with knots x and df
  ## this function is for x having unique values
  n = length(x);
  A = matrix(0, n, n);
  for(i in 1:n){
    y = rep(0, n); y[i]=1;
    yi = smooth.spline(x, y, df=df)$y;
    A[,i]= yi;
  }
  return(A)
}



S3 = smooth.matrix(x, df=3)
tmp=ns(x, df=3, intercept=TRUE)
H3 = tmp%*%solve(t(tmp)%*%tmp)%*%t(tmp);

## get the eigen value and eigen vector of the 
## smoother/projection matrix
eigen.S3 = eigen(S3)
eigen.H3 = eigen(H3)
v3 = eigen.S3$ve
dim(v3)

xx <- seq(-1,11,length=199)
fit.ss3<-smooth.spline(x, y,df = 3)
pss3 <- as.matrix(predict(fit.ss3, data.frame(x=xx))$y)
pss3
#4) B splines;knots at 0, -1 and 1 

bsplineBasis <-function (x, degree, innerknots, lowknot = min(x,innerknots), 
    highknot = max(x,innerknots)) {
    innerknots <- unique(sort (innerknots))
    knots <-c(rep(lowknot, degree + 1), innerknots, rep(highknot, degree + 1))
    n <- length (x)
    m <- length (innerknots) + 2 * (degree + 1)
    nf <- length (innerknots) + degree+1 
    basis <- rep (0,  n * nf)
    res <- splinebasis(d = degree,
      n = n, m =m, x = x, knots = knots, basis = basis)
    basis <- matrix (res$basis, n, nf)
    basis <- basis[,which(colSums(basis) > 0)]
    return (basis)
   
}


splinebasis<-function(d,n, m,x,knots,basis){
  mm = m
  dd = d
  nn = n
  k = mm - dd - 1
  for (i in 1:nn){
    ir = i + 1
    if (x[i] == knots[mm - 1]){
      basis[mindex(ir, k, nn) - 1] =  1
      for (j in 1:(k-1)){
        jr = j + 1
        basis[mindex (ir, jr, nn) - 1] = 0
      }
    } else {
      for (j in 1:k){
        jr = j + 1;
        basis [mindex (ir, jr, nn) - 1] = Mybs(mm, jr, dd + 1, x[i], knots);
      }
    }
  }
  return(list(basis=basis))
}

mindex<-function(i,j,nrow) {
  return ((j - 1) * nrow + i)
}

Mybs<-function(nknots,nspline,updegree,x,knots){
  if (updegree == 1) {
    if ((x >= knots[nspline - 1]) && (x < knots[nspline]))
      {y = 1}
    else {y = 0}
  }
  else {
    temp1 = 0
    if ((knots[nspline + updegree - 2] - knots[nspline - 1]) > 0)
    {temp1 = (x - knots[nspline - 1]) / (knots[nspline + updegree - 2] - knots[nspline - 1])}
    temp2 = 0
    if ((knots[nspline + updegree - 1] - knots[nspline]) > 0)
    {temp2 = (knots[nspline + updegree - 1] - x) / (knots[nspline + updegree - 1] - knots[nspline])}
    y1 = Mybs(nknots, nspline, updegree - 1, x, knots)
    y2 = Mybs(nknots, nspline + 1, updegree - 1, x, knots)
    y =  temp1 * y1 + temp2 * y2
  }
  return(y)
}
x<-rnorm(100,0,1)
b<-bsplineBasis (x, 3, c(-1,0,1))
fit <- lm( y ~ b)
print( summary(fit))
X<-cbind(rep(1,100),b)
Beta<-solve((t(X)%*%X))%*%t(X)%*%y
dim(Beta)
rowMeans(Beta)
#install.packages("splines2")
#### comparing with package ####
library(splines2)
bsMat <- bSpline(x, knots = c(-1,0,1), degree = 3, intercept = TRUE)
fit <- lm( y ~ bsMat)
print( summary(fit))
