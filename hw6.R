############# 5.2.1 ###############

# Sample from g
## function to generate a sample from g
g_sample <- function(n,theta) {
  c1 <- 2*gamma(theta)/(2*gamma(theta)+gamma(theta+0.5))
  uni <- runif(n)
  samples <- data.frame(spl=matrix(0,n,1))
  samples[which(uni<c1),1] <- rgamma(length(which(uni<c1)),shape = theta,scale = 1)
  samples[which(uni>=c1),1] <- rgamma(length(which(uni>=c1)),shape = theta+0.5,scale = 1)
  samples
}
## generate a sample with theta = 2
g_sample1 <- g_sample(10000,2)

## ture density
trued_g <- function (x,theta) {
  c1 <- 2*gamma(theta)/(2*gamma(theta)+gamma(theta+0.5))
  trued_g <- c1*dgamma(x,shape=theta,scale=1)+(1-c1)*dgamma(x,shape=theta+0.5,scale=1)
  trued_g
}
## plot
library(ggplot2)
ggplot(g_sample1) + stat_function(aes(col="true"), fun=trued_g, args = list(theta=2)) +
  geom_density(aes(spl, col='kernel')) + scale_color_manual(values = c("red", "blue")) +
  labs(title = expression(paste("kernel density vs. true density with ", theta, "=2")),
       x = "x", y = "g(x)") + theme(plot.title = element_text(hjust = 0.5))

# sample from f
## generate a sample from f
f_sample <- function(n,theta) {
  g_sample <- g_sample(n,theta)
  u_sample <- runif(n)
  c <- 1/(2*gamma(theta)+gamma(theta+0.5))
  f_sample <- matrix(0,n,1)
  for (i in 1:n) {
    qx <- sqrt(4 + g_sample[i,1])*g_sample[i,1]^(theta-1)*exp(-g_sample[i,1])
    gx <- (2*g_sample[i,1]^(theta-1)+g_sample[i,1]^(theta-0.5))*exp(-g_sample[i,1])
    if (u_sample[i] > qx*c/gx) next
      else {f_sample[i,1]<-g_sample[i,1]}
  }
  f_sample<-as.data.frame(f_sample[-which(f_sample==0),1])
  colnames(f_sample)<-"f"
  f_sample
}

f_sample1<-f_sample(10000,2)

qx<- function(x, theta) {
  qx<-sqrt(4 + x)*x^(theta-1)*exp(-x)
  qx
}
trued_f <- function (x,theta) {
  intgrl<-integrate(function(x) qx(x,theta),0,Inf)
  c_prime <-1/intgrl[["value"]]
  trued_f <- c_prime*sqrt(4 + x)*x^(theta-1)*exp(-x)
  trued_f
}

## plot
library(ggplot2)
ggplot() + stat_function(aes(col="true"), fun=trued_f, args = list(theta=2)) +
  geom_density(data=f_sample1,aes(f,col="kernel")) +
  scale_color_manual(values = c("red","blue")) +
  labs(title = expression(paste("kernel density vs. true density with ", theta, "=2")),
       x = "x", y = "f(x)") + theme(plot.title = element_text(hjust = 0.5))

############# 5.2.2 ###############
# Mixture Beta
# Sample from g
## function to generate a sample from g
g_sample <- function(n,theta,beta) {
  p1 <- beta/(beta+sqrt(3)*theta)
  uni <- runif(n)
  samples <- data.frame(spl=matrix(0,n,1))
  samples[which(uni<p1),1] <- rbeta(length(which(uni<p1)),shape1 = theta, shape2 = 1)
  samples[which(uni>=p1),1] <- rbeta(length(which(uni>=p1)),shape1 = 1, shape2 = beta)
  samples
}

# sample from f
## generate a sample from f
f_sample <- function(n,theta,beta) {
  g_sample <- g_sample(n,theta,beta)
  u_sample <- runif(n)
  f_sample <- matrix(0,n,1)
  for (i in 1:n) {
    qx <- g_sample[i,1]^(theta-1)/(1+g_sample[i,1]^2)+
      sqrt(2+g_sample[i,1]^2)*(1-g_sample[i,1])^(beta-1)
    hx <- g_sample[i,1]^(theta-1)+3^0.5*(1-g_sample[i,1])^(beta-1)
    if (u_sample[i] > qx/hx) next
    else {f_sample[i,1]<-g_sample[i,1]}
  }
  f_sample<-as.data.frame(f_sample[-which(f_sample==0),1])
  colnames(f_sample)<-"f"
  f_sample
}

f_sample1<-f_sample(10000,2,4)

qx<- function(x, theta, beta) {
  qx<-x^(theta-1)/(1+x^2)+sqrt(2+x^2)*(1-x)^(beta-1)
  qx
}
trued_f <- function (x,theta,beta) {
  intgrl<-integrate(function(x) qx(x,theta,beta),0,1)
  c_prime <-1/intgrl[["value"]]
  trued_f <- c_prime*(x^(theta-1)/(1+x^2)+sqrt(2+x^2)*(1-x)^(beta-1))
  trued_f
}

## plot
library(ggplot2)
ggplot() + stat_function(aes(col="true"), fun=trued_f, args = list(theta=2, beta=4)) +
  geom_density(data=f_sample1,aes(f,col="kernel")) +
  scale_color_manual(values = c("red","blue")) +
  labs(title = expression(paste("kernel density vs. true density with ", theta, 
                                "=2 and ", beta, "=4")),
       x = "x", y = "f(x)") + theme(plot.title = element_text(hjust = 0.5))

# Separate Beta

# Sample from g
## function to generate a sample from g
g_sample1 <- function(n,theta,beta) {
  g_sample1 <- data.frame(rbeta(n,shape1=theta,shape2=1))
  g_sample1
}

g_sample2 <- function(n,theta,beta) {
  g_sample2 <- data.frame(rbeta(n,shape1=1,shape2=beta))
  g_sample2
}

# sample from f
## generate a sample from f
f_sample_sep <- function(n,theta,beta) {
  u_sample <- runif(n)
  n1 <- length(which(u_sample<0.5))
  g_sample1 <- g_sample1(n1,theta,beta)
  u_sample1 <- runif(n1)
  f_sample1 <- matrix(0,n1,1)
  for (j in 1:n1) {
    qx1 <- g_sample1[j,1]^(theta-1)/(1+g_sample1[j,1]^2)
    hx1 <- g_sample1[j,1]^(theta-1)
    if (u_sample1[j] > qx1/hx1) next
     else {f_sample1[j,1]<-g_sample1[j,1]}
   }
  n2 <- length(which(u_sample>=0.5))
  g_sample2 <- g_sample2(n2,theta,beta)
  u_sample2 <- runif(n2)
  f_sample2 <- matrix(0,n2,1)
  for (k in 1:n2) {
    hx2 <- 3^0.5*(1-g_sample2[k,1])^(beta-1)
    qx2 <- sqrt(2+g_sample2[k,1]^2)*(1-g_sample2[k,1])^(beta-1)
    if (u_sample2[k] > qx2/hx2) next
    else {f_sample2[k,1]<-g_sample2[k,1]}
    }
  f_sample<-rbind(f_sample1,f_sample2)
  f_sample<-as.data.frame(f_sample[-which(f_sample==0),1])
  colnames(f_sample)<-"f"
  f_sample
}  

f_sample_sep<-f_sample_sep(10000,2,4)

## plot
library(ggplot2)
ggplot() + stat_function(aes(col="true"), fun=trued_f, args = list(theta=2, beta=4)) +
  geom_density(data=f_sample_sep,aes(f,col="kernel")) +
  scale_color_manual(values = c("red","blue")) +
  labs(title = expression(paste("kernel density vs. true density with ", theta, 
                                "=2 and ", beta, "=4")),
       x = "x", y = "f(x)") + theme(plot.title = element_text(hjust = 0.5))


