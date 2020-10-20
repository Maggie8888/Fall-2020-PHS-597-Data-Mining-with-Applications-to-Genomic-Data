
rm(list=ls())
library(MASS)
library(mvtnorm)
library(ggplot2)
Mysvm <- function(X, Y, C, kernelFunction,tol, max.iter, sigma=0.5,d=4){
   
  ## data parameter
  m <- nrow(X)
 
  ## labels
  ## Map 0 to -1
  Y[Y==0] <- -1
  
  # variables
  alphas <- rep(0, m)
  b <- 0
  E <- rep(0, m)
  iter <- eta <- L <- H <- 0
  
  ## transformations of original data to map into new space
  if(kernelFunction == "linearKernel") {
    K <- X %*% t(X)
  } else if (kernelFunction == "gaussianKernel") {
    K <- matrix(0, ncol=m, nrow=m)
    for (i in 1:m) {
      for (j in i:m) {
        K[i,j] <- gaussianKernel(X[i,], X[j,],sigma)
        K[j,i] <- K[i,j]
      }
    }
  }
  
  else if (kernelFunction == "polynomialKernel") {
    K <- matrix(0, ncol=m, nrow=m)
    for (i in 1:m) {
      for (j in i:m) {
        K[i,j] <- polynomialKernel(X[i,], X[j,],d)
        K[j,i] <- K[i,j]
      }
    }
  }
  
  while (iter < max.iter) {
    num_changed_alphas <- 0
    for (i in 1:m){
      E[i] <- b + sum(alphas * Y * K[,i]) - Y[i]
      if( (Y[i]*E[i] < -tol & alphas[i] < C) || (Y[i] >  tol & alphas[i] > 0) ) {
        ## if the error E[i] is large
        ## the alpha corresponding to this data instance can be optimized
        
        j <- ceiling(m * runif(1))
        while (j == i) { # make sure i != j
          j <- ceiling(m * runif(1))
        }
        
        E[j] <- b + sum(alphas * Y * K[,j]) - Y[j]
        
        ## save old alphas
        alpha.i.old <- alphas[i]
        alpha.j.old <- alphas[j]
        
        if (Y[i] == Y[j]) {
          L <- max(0, alphas[j] + alphas[i] -C)
          H <- min(C, alphas[j] + alphas[i])
        } else {
          L <- max(0, alphas[j] - alphas[i])
          H <- min(C, C + alphas[j] - alphas[i])
        }
        
        if (L == H) {
          ## continue to next i
          next
        }
        
        ## compute eta
        eta <- 2 * K[i,j] - K[i,i] - K[j,j]
        if (eta >= 0) {
          ## continue to next i
          next
        }
        
        ## compute and clip new value for alpha j
        alphas[j] <- alphas[j] - Y[j] * (E[i] - E[j])/eta
        ## clip
        alphas[j] = min(H, alphas[j])
        alphas[j] = max(L, alphas[j])
        
        ## check if change in alpha is significant
        if(abs(alphas[j] - alpha.j.old) < tol) {
          alphas[j] <- alpha.j.old
          next
        }
        
        ## determine value for alpha i
        alphas[i] <- alphas[i] + Y[i] * Y[j] * (alpha.j.old - alphas[j])
        
        ## compute b1 and b2
        b1 <- b - E[i] +
          - Y[i] * (alphas[i] - alpha.i.old) * K[i,j] +
          - Y[j] * (alphas[j] - alpha.j.old) * K[i,j]
        b2 <- b - E[j] +
          - Y[i] * (alphas[i] - alpha.i.old) * K[i,j] +
          - Y[j] * (alphas[j] - alpha.j.old) * K[i,j]
        
        ## compute b
        if ( alphas[i] > 0 & alphas[i] < C) {
          b <- b1
        } else if (alphas[j] > 0 & alphas[j] < C) {
          b <- b2
        } else {
          b <- (b1+b2)/2
        }
        
        num_changed_alphas <- num_changed_alphas + 1
      }
    }
    
    if (num_changed_alphas == 0) {
      iter <- iter + 1
    } else {
      iter <-  0
    }
    
  }
  
  
  idx <- alphas > 0
  return(list(
      X=X[idx,],
      y=Y[idx],
      kernelFunction=kernelFunction,
      b=b,
      alphas=alphas[idx],
      w=t(alphas * Y) %*% X
  ))
}


gaussianKernel <- function(x1,x2, sigma) {
  sim <- exp(-sum((x1-x2)^2)/(2*sigma^2))
  return(sim)
}



polynomialKernel <- function(x1,x2, d) {
  sim <-(t(x1)%*%x2+1)^d
  return(sim)
}



svmPredict <- function(model, X,kernelFunction,sigma,d) {
  ## returns a vector of predictions using a trained SVM model
  ## X is a m x n matrix where there each example is a row.
  ## model is a svm model returned from svm.
  m <- nrow(X)
  p <- rep(0, m)
  pred <- rep(0, m)
  w<-model$w
  b<-model$b
  if (kernelFunction == "linearKernel") {
    p <- X %*% t(w) +b
  } else if (kernelFunction == "gaussianKernel") {
    
    alphas <- model$alphas
    for (i in 1:m) {
      prediction <- 0
      for (j in 1:nrow(model$X)) {
        prediction <- prediction +
          alphas[j] * model$y[j] *
          gaussianKernel(X[i,], model$X[j,],sigma)
      }
      p[i] <- prediction + model$b
    }
    
  }
  else if (kernelFunction == "polynomialKernel") {
    
    alphas <- model$alphas
    for (i in 1:m) {
      prediction <- 0
      for (j in 1:nrow(model$X)) {
        prediction <- prediction +
          alphas[j] * model$y[j] *
          polynomialKernel(X[i,], model$X[j,],d)
      }
      p[i] <- prediction + model$b
    }
    
  }
  
  
  pred[p >= 0] <- 1
  pred[p < 0] <- -1
  return(pred)
}


## equivalent with meshgrid in octave
meshgrid <- function(a, b) {
  list(
    x <- outer(b*0, a, FUN="+"),
    y <- outer(b, a*0, FUN="+")
  )
}

Myplot<-function(x, X, y, type, title="", xlab="", ylab="",sigma,d){
  model <- x
  match.arg(type, c("linearKernel", "gaussianKernel","polynomialKernel"))
  plot_data<-as.data.frame(X)
  plot_data<-cbind(plot_data,as.factor(y))
  colnames(plot_data)<-c("X1","X2","y")
  p <- ggplot(data=plot_data, aes_string(x=colnames(plot_data)[1], y=colnames(plot_data)[2])) +
    geom_point(aes(colour=y))
 
  w <- model$w
  b <- model$b
  if(type == "linearKernel") {
    p <- p+geom_abline(intercept = -b/w[2],
                       slope = -w[1]/w[2],
                       colour = "red")
  } else {
    xr <- seq(min(X[,1]), max(X[,1]), length.out=200)
    yr <- seq(min(X[,2]), max(X[,2]), length.out=200)
    mg <- meshgrid(xr, yr)
    X1 <- mg[[1]]
    X2 <- mg[[2]]
    vm<-as.data.frame(cbind(X1,X2))
    vals <- matrix(0, ncol=ncol(X1), nrow=nrow(X1))
    z<-vector()
    for (i in 1:ncol(X1)){
      thisX <- cbind(X1[,i], X2[,i])
      vals[,i] <-  svmPredict(model, thisX,kernelFunction=type,sigma,d)
    }
    
    vm <- data.frame(x=as.vector(X1),
                     y= as.vector(X2),
                    z=as.vector(vals))
   
    #p<-ggplot(vm, aes(x =x, y = y, shape = factor(z))) +
    #  geom_point(colour = "grey90", size = 1.5) + 
    # geom_density_2d(aes(colour = z))+stat_contour(breaks = c(0))
    
    p<-ggplot(vm, aes(x =x, y = y, z=z),shape = factor(z)) +
      geom_point(aes(colour=z)) + 
    stat_contour(breaks = c(0))
    
    
  }
  
  p <- p + xlab(xlab) + ylab(ylab) + ggtitle(title)
  
  print(p)
}
set.seed(1234567)
#### linear seperable cases ####
x1=rmvnorm(100,mean=c(1.5,1.5),sigma=0.01*diag(2))
x2=rmvnorm(100,mean=c(2,2),sigma=0.001*diag(2))


X<-rbind(x1,x2)
colnames(X)<-c("1","2")
y=matrix(0,nrow=200,ncol=1)
y[1:100,1]=1
y[101:200,1]=-1
C=0.1

model<-Mysvm(X, y, C,"linearKernel",
                tol=1e-4, max.iter=20,sigma=0.5,d=4)
pred<-svmPredict(model,X,"linearKernel",sigma=0.5,4)
Myplot(model, X, y,"linearKernel", title="", xlab="", ylab="",sigma=0.5,d=4)

model<-Mysvm(X, y, C,"gaussianKernel",
                tol=1e-4, max.iter=20,sigma=0.5,d=4)
pred<-svmPredict(model,X,"gaussianKernel",sigma=0.5,4)
Myplot(model, X, y,"gaussianKernel", title="", xlab="", ylab="",sigma=0.5,d=4)

model<-Mysvm(X, y, C,"polynomialKernel",
                     tol=1e-4, max.iter=20, sigma=0.5,d=4)
pred<-svmPredict(model,X,"polynomialKernel",sigma=0.5,4)
Myplot(model, X, y,"polynomialKernel", title="", xlab="", ylab="",sigma=0.5,d=4)
##################### using existed function #################
library(e1071)

m <- svm(y ~ X, scale = FALSE, kernel = "linear")
coef(m)
plot_data<-as.data.frame(X)
plot_data<-cbind(plot_data,as.factor(y))
colnames(plot_data)<-c("X1","X2","y")
p <- ggplot(data=plot_data, aes_string(x=colnames(plot_data)[1], y=colnames(plot_data)[2])) +
  geom_point(aes(colour=y))+geom_abline(intercept = -coef(m)[1]/coef(m)[3],
                                        slope = -coef(m)[2]/coef(m)[3],
                                        colour = "red")

p

##### linear non-seperable case ####
x1=rmvnorm(100,mean=c(1,1),sigma=0.5*diag(2))
x2=rmvnorm(100,mean=c(0,0),sigma=0.01*diag(2))
y=matrix(0,nrow=200,ncol=1)

X<-rbind(x1,x2)
colnames(X)<-c("1","2")

y<-ifelse(X[,1]>0|X[,2]>0,1,-1)

C=100

model<-Mysvm(X, y, C,"linearKernel",
           tol=1e-4, max.iter=20,sigma=0.5,d=4)
pred<-svmPredict(model,X,"linearKernel",sigma=0.5,4)
Myplot(model, X, y,"linearKernel", title="", xlab="", ylab="",sigma=0.5,d=4)

model<-Mysvm(X, y, C,"gaussianKernel",
           tol=1e-4, max.iter=20,sigma=0.5,d=4)
pred<-svmPredict(model,X,"gaussianKernel",sigma=0.5,4)
Myplot(model, X, y,"gaussianKernel", title="", xlab="", ylab="",sigma=0.5,d=4)

model<-Mysvm(X, y, C,"polynomialKernel",
           tol=1e-4, max.iter=20, sigma=0.5,d=4)
pred<-svmPredict(model,X,"polynomialKernel",sigma=0.5,4)
Myplot(model, X, y,"polynomialKernel", title="", xlab="", ylab="",sigma=0.5,d=4)


############## using existed function ########

m <- svm(y ~ X, scale = FALSE, kernel = "radial",gamma=0.5,degree=4,cost=1)
m$decision.values
predict(m,X)

xr <- seq(min(X[,1]), max(X[,1]), length.out=200)
yr <- seq(min(X[,2]), max(X[,2]), length.out=200)
mg <- meshgrid(xr, yr)
X1 <- mg[[1]]
X2 <- mg[[2]]
vm<-as.data.frame(cbind(X1,X2))
vals <- matrix(0, ncol=ncol(X1), nrow=nrow(X1))
z<-vector()
for (i in 1:ncol(X1)){
  thisX <- cbind(X1[,i], X2[,i])
  vals[,i] <-  predict(m, thisX)
}

vm <- data.frame(x=as.vector(X1),
                 y= as.vector(X2),
                 z=as.vector(vals))


p<-ggplot(vm, aes(x =x, y = y, z=z),shape = factor(z)) +
  geom_point(aes(colour=z)) + 
  stat_contour(breaks = c(0))
p

