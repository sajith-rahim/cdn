#---------------Basic Sampling Algorithm-------------------
#-----------Probability Integral Transformation------------

# return n samples based on the inverse of the target cdf
inv.transform <- function(inv.f, n) {
  
# Non-vectorized version (for explanation purposes only)  
#
#   result.sample <- rep(NA,n)
#   
#   for (i in 1:n) {
#     u <- runif(1,0,1)              # step 1
#     result.sample[i] <- inv.f(u)   # step 2
#   }
#   
#   result.sample

# Vectorized version
  inv.f(runif(n,0,1))
}


#EG:1 fx(x) = 3x^2
#FX(x) = x^3
#FX^-1(u) = u^(1/3)


inv.f <- function(u) u^(1/3)

vals <- inv.transform(inv.f, 5e4)

# Plotting
hist(vals, breaks=50, freq=FALSE, main=expression("Sample vs true Density [ f(x)=" ~3*x^2~"]"))
curve(3*x^2, 0, 1, col="red", lwd=2)




#EG:1 fx(x) = 1/3(x^2)
#FX(x) = 1/9(x^3+1)
#FX^-1(u) = (9u-1)^1/3


#   inv.f <- function(u) (9*u-1)^(1/3)
inv.f <- function(u) ifelse((9*u-1)>=0,  (9*u-1)^(1/3),  -(1-9*u)^(1/3))

vals <- inv.transform(inv.f, 5e4)

# Plotting
hist(vals, breaks=50, freq=FALSE, main="Sample vs true Density")
curve((1/3)*x^2, -1, 2, col="red", lwd=2, add=T)