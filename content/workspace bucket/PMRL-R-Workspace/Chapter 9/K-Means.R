#PMRL Chapter 9

#----------------Mixture Models------------
# Refer:
#    Gaussian Mixture Models
#    K-Means 9.1
#    
#----------------------------------------------------------
#--------K-Means-------------------------------------------


install.packages("MASS")
library(MASS)

# generating n datapoints from a mixture of K d dim Gaussians
# k  : the respective datapoint classes
# mu : kxd matrix with means
# sig: kxdxd matrix with dxd covariate matrices
gen.mix <- function(n, k, mu, sig) {
  
  d <- length(mu[1,])  # number of dimensions
  result <- matrix(rep(NA,n*d), ncol=d)
  colnames(result) <- paste0("X",1:d)
  
  for(i in 1:n) {
    result[i,] <- mvrnorm(1, mu = mu[k[i],], Sigma=sig[,,k[i]])
  }
  
  result
}

set.seed(101)
n <- 360

mu <- matrix(c(14.0,4.0,
               15.0,5.0,
               16.5,5.0), ncol=2, byrow=T)

sigs <- array(rep(NA,2*2*3), c(2,2,3))  # 3D matrix
sigs[,,1] <- matrix(c(.25, .21, .21,.25), nrow=2, byrow=TRUE)
sigs[,,2] <- matrix(c(.25,-.21,-.21,.25), nrow=2, byrow=TRUE)
sigs[,,3] <- matrix(c(.25, .21, .21,.25), nrow=2, byrow=TRUE)

pi      <- c(.2,.5,.3)                           # mixing coeffs
classes <- sample(1:3, n, replace=TRUE, prob=pi)

d <- gen.mix(n, classes, mu, sigs)
plot(d, col=c("red","green","blue")[classes], xlab="X1", ylab="X2", pch=19)



k_means <- function(dataset, K, max_iter=100) {
  
  # get the dataset classification given the current indicators
  get_classes <- function(rnk)
    apply(rnk,1,function(row) which.max(row))
  
  d      <- ncol(dataset)                                # number of dimensions
  N      <- nrow(dataset)                                # number of samples
  ranges <- sapply(1:d, function (i) range(dataset[,i])) # the ranges for each dimension
  
  # generate K initial random cluster centers (each center is a row vector)
  mu <- t(replicate(K,sapply(1:d, function(i) runif(1,ranges[1,i], ranges[2,i]))))
  
  # indicators (each row consists of 0...1...0, ie, it's a 1-of-K coding scheme)
  rnk <- matrix(rep(0,K*n), ncol=K)
  old_classes <- get_classes(rnk)
  
  for(it in 1:max_iter) {
    
    # update indicators for each datapoint
    for(n in 1:N) {
      distances <- sapply(1:K, function(k) norm(as.matrix(dataset[n,]-mu[k,]),"F"))
      rnk[n,]   <- rep(0,K)
      rnk[n,which.min(distances)] <- 1
    }
    
    classes <- get_classes(rnk)
    if (all(old_classes == classes)) # convergence achieved?
      break
    else 
      old_classes <- classes
    
    # update centers given the updated indicators
    for(k in 1:K) {
      mu[k,]  <- rnk[,k] %*% dataset / sum(rnk[,k])
    }
  }
  
  list(mu=mu, pred=classes)
}


set.seed(101)
result <- k_means(d,3)  # set clustering to 3 classes

plot(d, col=c("red","green","blue")[result$pred], xlab="X1", ylab="X2", pch=19)
points(result$mu, pch=3, lwd=5)  # plot the centers
