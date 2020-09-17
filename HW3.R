
## illustration of a scenario where linear regression on the indicator matrix cannot distinguish classes 
rm(list=ls())
library(MASS)
library(mvtnorm)
library(ggplot2)
set.seed(12345678)
x1=rmvnorm(100,mean=c(1,1),sigma=0.01*diag(2))
x2=rmvnorm(100,mean=c(2,2),sigma=0.1*diag(2))
x3=rmvnorm(100,mean=c(3,3),sigma=0.01*diag(2))

plot_data<- rbind(x1,x2,x3)
plot_data<-as.data.frame(plot_data)
plot_data<-cbind(plot_data,c(rep(c("A","B","C"),each=100)))
colnames(plot_data)<-c("x","y","group")

p<-ggplot(plot_data,aes(x,y, shape = factor(group))) +  
  geom_point(aes(colour = factor(group)), size = 4) +
  geom_point(colour = "grey90", size = 1.5)


x=cbind(rep(1,300),rbind(x1,x2,x3));
y=matrix(0,nrow=300,ncol=3);
y[1:100,1]=1
y[101:200,2]=1;
y[201:300,3]=1;
#obtain regression coefficient;
b=ginv(t(x)%*%x)%*%(t(x)%*%y);
#predict values of Y;
yhat=x%*%b;
#The classification is done using the predicted values of yhat;
y.pred=apply(yhat,1,which.max);
# draw decision boundary using the estimated coefficient of b; 
b12=b[,1]-b[,2];
b23=b[,2]-b[,3];
b13=b[,1]-b[,3];
# draw decision boundary on top of the graph; 
p + geom_abline(intercept = -b12[1]/b12[3], slope = -b12[2]/b12[3], color="red", 
            linetype="dashed", size=1)+
  geom_abline(intercept = -b23[1]/b23[3], slope = -b23[2]/b23[3], color="blue", 
              linetype="dotted", size=1)+
  geom_abline(intercept = -b13[1]/b13[3], slope = -b13[2]/b13[3], color="green", 
             linetype="dotted", size=1)

##############################################################################
#LDA illustration:
#estimate the center for each class;
x1=rmvnorm(100,mean=c(1,1),sigma=0.01*diag(2))
x2=rmvnorm(100,mean=c(2,2),sigma=0.1*diag(2))
x3=rmvnorm(100,mean=c(3,3),sigma=0.01*diag(2))

x=rbind(x1,x2,x3);n=nrow(x)
mu1=colMeans(x1);mu2=colMeans(x2);mu3=colMeans(x3);
pi1=nrow(x1)/(nrow(x1)+nrow(x2)+nrow(x3));
pi2=nrow(x2)/(nrow(x1)+nrow(x2)+nrow(x3));
pi3=nrow(x3)/(nrow(x1)+nrow(x2)+nrow(x3));
sigma=1/(n-3)*((t(x1)-mu1)%*%(x1-mu1)+(t(x2)-mu2)%*%(x2-mu2)+(t(x3)-mu3)%*%(x3-mu3))
#decision boundary:
b12=ginv(sigma)%*%(mu2-mu1)
a12=-1/2*(mu1+mu2)%*%ginv(sigma)%*%(mu2-mu1)


b13=ginv(sigma)%*%(mu3-mu1)
a13=-1/2*(mu1+mu3)%*%ginv(sigma)%*%(mu3-mu1)

b23=ginv(sigma)%*%(mu3-mu2)
a23=-1/2*(mu2+mu3)%*%ginv(sigma)%*%(mu3-mu2)

p + geom_abline(intercept = -a12/b12[2], slope = -b12[1]/b12[2], color="red", 
                linetype="dashed", size=1.5)+
  geom_abline(intercept = -a13/b13[2], slope = -b13[1]/b13[2], color="blue", 
              linetype="dotted", size=1.5)+
  geom_abline(intercept = -a23/b23[2], slope = -b23[1]/b23[2], color="green", 
              linetype="dotted", size=1.5)


