#PMRL Chapter 9

#----------------Mixture Models EM------------
# Refer:
#    Gaussian Mixture Models
#    EM 9.2.2
#    
#----------------------------------------------------------
#--------GMM EM-------------------------------------------


install.packages("MASS")
library(MASS)
install.packages('mvtnorm')
library(mvtnorm)

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



plot(d, col="black", xlab="X1", ylab="X2", pch=19)



#library(mvtnorm)

em_gaussian_mix <- function(dataset, K, max_iter=100, epsilon=1e-3) {
  
  # get the dataset classification given the current indicators
  get_classes <- function(gammak)
    apply(gammak,1,function(row) which.max(row))
  
  d      <- ncol(dataset)                                # number of dimensions
  N      <- nrow(dataset)                                # number of samples
  ranges <- sapply(1:d, function (i) range(dataset[,i])) # the ranges for each dimension
  
  # initial values
  pik <- rep(1/K,K)
  muk <- t(replicate(K,sapply(1:d, function(i) runif(1,ranges[1,i], ranges[2,i]))))
  Sigmas <- array(rep(NA,2*2*3), c(2,2,3)) 
  for (k in 1:K)
    Sigmas[,,k] <- diag(d)
  gammak <- matrix(rep(0,K*N),ncol=K) # the responsabilities
  old_gammak <- gammak
  
  # EM steps
  for(it in 1:max_iter) {
    
    # Expectation step: compute responsibilities
    
    for (k in 1:K) {
      gammak[,k] <- apply(dataset, 1, 
                          function(xi) {
                            pik[k] * dmvnorm(xi,muk[k,], Sigmas[,,k])
                          })
    }
    gammak <- t(apply(gammak, 1, function(row) row/sum(row)))
    
    if (sum(abs(gammak - old_gammak)) < epsilon) # converged?
      break
    else 
      old_gammak <- gammak
    
    # Maximization step: maximize the expected value wrt parameters theta
    
    Nk  <- sapply(1:K, function (k) sum(gammak[,k]))
    pik <- Nk/N
    for (k in 1:K) {
      muk[k,]     <- apply(gammak[,k] * dataset,2,sum) / Nk[k]
      Sigmas[,,k] <- diag(d) * 0 # reset
      for(n in 1:N) {
        Sigmas[,,k] <- Sigmas[,,k] + 
          gammak[n,k]* (dataset[n,]-muk[k,])%*%t(dataset[n,]-muk[k,])
      }
      Sigmas[,,k] <- Sigmas[,,k] / Nk[k]  
    }
  }
  
  list(mu=mu, Sigmas=Sigmas, gammak=gammak, pred=get_classes(gammak))
}

set.seed(101)
result <- em_gaussian_mix(d,3)  # set clustering to 3 classes

plot(d, col=c("red","green","blue")[result$pred], xlab="X1", ylab="X2", pch=19)

#RGB to represent mixture
plot(d, col=rgb(result$gammak[,1], result$gammak[,2], result$gammak[,3]), xlab="X1", ylab="X2", pch=19)