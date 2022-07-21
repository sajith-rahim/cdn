#PMRL Chapter 4

#----------------Descriminant Functions-----------------------
# Refer:
#   Descriminant Functions
#   2Class Descriminanant and LS for Classification
# Issue:
#   Poor Result
#----------------------------------------------------------
#----------------Descriminant Functions--------------------


#----------------Demo Data Setup--------------------

library(MASS)


genetrate.mixture <- function(n, classes, mu, sig) {
  # generating n datapoints from a mixture of K Gaussians with dimensions d
  # classes  : the respective datapoint classes
  # mu : kxd matrix with means
  # sig: kxdxd matrix with dxd covariate matrices
  
  d <- length(mu[1,])  # number of dimensions
  result <- matrix(rep(NA,n*d), ncol=d)
  colnames(result) <- paste0("X",1:d)
  
  #mvrnorm From MASS - Simulate From A Multivariate Normal Distribution 
  #Produces one or more samples from the specified multivariate normal distribution.
  for(i in 1:n) {
    result[i,] <- mvrnorm(1, mu = mu[classes[i],], Sigma=sig[,,classes[i]])
  }
  
  result
}


set.seed(101)
n <- 100

mu <- matrix(c(4.0,4.0,6.5,  5), ncol=2, byrow=T)
#       [,1] [,2]
# [1,]  4.0    4
# [2,]  6.5    5

sigs <- array(rep(NA,2*2*2), c(2,2,2))  # 3D matrix
sigs[,,1] <- matrix(c(.25, .21, .21,.25), nrow=2, byrow=TRUE)
sigs[,,2] <- matrix(c(.25, .21, .21,.25), nrow=2, byrow=TRUE)
# , , 1
# 
# [,1] [,2]
# [1,] 0.25 0.21
# [2,] 0.21 0.25
# 
# , , 2
# 
# [,1] [,2]
# [1,] 0.25 0.21
# [2,] 0.21 0.25

pi <- c(.6,.4) # mixing coefficientss
classes <- sample(1:2, n, replace=TRUE, prob=pi)

mydata <- genetrate.mixture(n, classes, mu, sigs)
mydata <- cbind(mydata, C=classes)
head(mydata)

plot(mydata, col=c("red","blue")[classes], xlab="X1", ylab="X2", pch=19)


#------Least Squares for Classification------

compute_W <- function(X, C) {
  solve(t(X) %*% X) %*% t(X) %*% C
}

X <- matrix(cbind(rep(1,length(mydata[,1])), mydata[,1:2]), ncol=3)
C <- matrix(mydata[,3], ncol=1)
W <- compute_W(X, C)


# given a W, create a discriminant function defined by it
discriminant_factory <- function(W) {
  function(x) t(W) %*% c(1,x)
}

f <- discriminant_factory(W)  # this problem's discriminant function


set.seed(231)
n <- 200
unknown_class <- sample(1:2, n, replace=TRUE, prob=pi)
newdata       <- genetrate.mixture(n, unknown_class, mu, sigs)
plot(newdata, pch=19)

#threshold is -w0 (check?)
prediction <- ifelse(apply(newdata, 1, f)<=1.5, 1, 2) # apply threshold
plot(newdata, col=c("red","blue")[prediction], xlab="X1", ylab="X2", pch=19)

# draw decision line: w0 + w1x1 + w2x2 =0 (this should be the line)
abline((-W[1])/W[3], -W[2]/W[3], col = "red" ,lty=4)  



#?
# draw decision line: w0 + w1x1 + w2x2 == 1.5
abline((1.5-W[1])/W[3], -W[2]/W[3], lty=3)  


#Contingency Table
table(unknown_class, prediction)
