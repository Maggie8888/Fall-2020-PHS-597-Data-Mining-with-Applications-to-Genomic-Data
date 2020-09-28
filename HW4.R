
rm(list=ls())

library(MASS)
library(mvtnorm)
library(ggplot2)
set.seed(12345678)
x1=rmvnorm(100,mean=c(1,1),sigma=0.01*diag(2))
x2=rmvnorm(100,mean=c(2,2),sigma=0.1*diag(2))
x3=rmvnorm(100,mean=c(3,3),sigma=0.01*diag(2))
X=rbind(x1,x2,x3)
y=c(rep(1,100),rep(2,100),rep(3,100))
s <- sort(sample(nrow(X), 4))
print(X[s, ])
print(y[s])
tol <- 1e-04
X <- as.matrix(X)
n <- nrow(X)
p <- ncol(X)
g <- as.factor(y)
lev <- lev1 <- levels(g)
counts <- as.vector(table(g))
proportions <- counts/n
ng <- length(proportions)
names(counts) <- lev1 # `names(prior)` excised
# group counts and proportions of the total:
print(cbind(counts, proportions))

## drop attributes to avoid e.g. matrix() methods
group.means <- tapply(c(X), list(rep(g, p), col(X)), mean)
f1 <- sqrt(diag(var(X - group.means[g,  ])))
# scale columns to unit variance before checking for collinearity
scaling <- diag(1/f1, p)

# group centroids
print(group.means)

# within-group standard deviations
print(f1)
# inverses of within-group standard deviations
scaling0 <- scaling
print(scaling0)
# within-group correlation matrix
print(cor(X - group.means[g, ]))
# within-group covariance matrix after standardization
print(cov((X - group.means[g, ]) %*% scaling))

fac <- 1/n # excise conditional case `method == "moment"`
XX <- sqrt(fac) * (X - group.means[g,  ]) %*% scaling
X.s <- svd(XX, nu = 0L)
rank <- sum(X.s$d > tol)
scaling <- scaling %*% X.s$v[,1L:rank] %*% diag(1/X.s$d[1L:rank],rank)

# standardized within-group centered data
Xhat <- (X - group.means[g,  ]) %*% scaling0
# eigendecomposition of covariance
E <- eigen(cov(Xhat))
# sphering matrix (alternative calculation)
scaling0 %*% E$vectors %*% diag(1/sqrt(E$values))
# sphering matrix
scaling1 <- scaling
print(scaling1)
# within-group covariance matrix after sphering
print(cov((X - group.means[g, ]) %*% scaling))
xbar <- colSums(proportions %*% group.means) # sub `proportions` for `prior`
fac <- 1/(ng - 1) # excise conditional case `method == "mle"`
XX <- sqrt((n * proportions)*fac) * # sub `proportions` for `prior`
  scale(group.means, center = xbar, scale = FALSE) %*% scaling
X.s <- svd(XX, nu = 0L)
rank <- sum(X.s$d > tol * X.s$d[1L])
scaling <- scaling %*% X.s$v[, 1L:rank]
# excise conditional case `is.null(dimnames(x))`
dimnames(scaling) <- list(colnames(X), paste("LD", 1L:rank, sep = ""))
dimnames(group.means)[[2L]] <- colnames(X)
# number of dimensions in the sphering-transformed group centroid SVD
print(rank)

# discriminant coefficients
scaling2 <- scaling
print(scaling2)
# discriminant coordinates of centered group centroids
print((group.means - matrix(1, ng, 1) %*% t(xbar)) %*% scaling2)
# as recovered using `predict.lda()`
lda.model <- lda(X, grouping = y)
print(predict(lda.model, as.data.frame(lda.model$means))$x)
predict(lda.model,X)$class

# center data around weighted means & transform
means <- colSums(lda.model$prior * lda.model$means)
mod <- scale(X, center = means, scale = FALSE) %*% 
  lda.model$scaling
# visualize the features in the two LDA dimensions
plot.df <- data.frame(mod, "Outcome" = y)
library(ggplot2)
ggplot(plot.df, aes(x = LD1, y = LD2, color = Outcome)) + geom_point()
################## QDA ######################
QuaDA<- function(Input,Response){
  Input=as.matrix(Input)
  N=nrow(Input)
  p=ncol(Input)
  K=length(unique(Response))
  classes=as.numeric(levels(as.factor(Response)))
  Mu_k=matrix(0,nrow = p,ncol = K)#centroid vectors
  Pi_k=vector(length= K)# classes a priori
  sigma_k=list()#data covariances matrix per class
  sigmaMinus1_k=list()
  det_k=vector(length = K)
  for (i in 1:K){
    sigma_k[[i]]=matrix(0,ncol=p,nrow = p)
    tmp=Input[Response==classes[i],]
    N0=nrow(tmp)
    Pi_k[i]=N0/N
    Mu_k[,i]=colMeans(tmp)
    for(j in 1:N0){
      sigma_k[[i]]=sigma_k[[i]]+(tmp[j,]-Mu_k[,i])%*%t(tmp[j,]-Mu_k[,i])
    }
    sigma_k[[i]]=(1/(N0-1))*sigma_k[[i]]
    sigmaMinus1_k[[i]]=solve(sigma_k[[i]])
    det_k[i]=det(sigma_k[[i]])
  }
  #discriminant matrix
  deltaTrain=matrix(0,ncol = K,nrow = N)
  for (i in 1:K){
    deltaTrain[,i]=apply(Input,1,function(x) -0.5*t(x-Mu_k[,i])%*%sigmaMinus1_k[[i]]%*%(x-Mu_k[,i]))
    deltaTrain[,i]=deltaTrain[,i]-0.5*log(det_k[i])+log(Pi_k[i])
  }
  #prediction
  yHatTrain=apply(deltaTrain,1,which.max)
  return(list(yhat=yHatTrain,sigmaMinus1_k=sigmaMinus1_k,Mu_k=Mu_k,det_k=det_k,Pi_k=Pi_k))
}

QuaDA(X,y)
