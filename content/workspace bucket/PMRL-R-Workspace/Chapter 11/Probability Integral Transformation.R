#-----------Probability Integral Transformation------------
# A continous random variable transforemd by it's own CDF will always
#follow a U(0,1) distribution
#We can thus transform a U(0,1) rv with inverse of the CDF
# to get the RV with that CDF
#    
#----------------------------------------------------------



n <- 1000
lambda <- 5
#x ~ f(x) = ?? {e}^{- ?? x} with ?? =5
x <- rexp(n, lambda)
#y <- cdf of x
y <- 1 - exp(-lambda*x)
F_y <- cumsum(table(y))/sum(table(y))
#F_y is the distribution function of U(0,1)
par(mfrow=c(1,3))
hist(x, col='blue')
plot(F_y, col='green', pch=19, xlab='y', ylab=expression('F'['Y']), main=expression(paste('Y=1-e'^ paste('-',lambda,'X'))))
plot(ecdf(y), col='red', pch=19, xlab='y', ylab='ECDF(y)')