#PMRL Chapter 3

#----------------Linear Basis Models-----------------------
# Refer:
#   Different Basis Functions and Locality - Polynomial Guassian Sigmoid etc
#   Maximum Likelihood t = y(x,w) + e ; e ~ Norm Derivation
#   Geometric interpretation of SSE with Phi
#----------------------------------------------------------
#----------------Linear Basis Models MLE-------------------
  

# compute weights for a given basis phi, using maximum likelihood
compute_w <- function(X, Y, phi) {
  Phi <- sapply(phi, function(basisFn) basisFn(X))  # make design matrix
  solve(t(Phi) %*% Phi) %*% t(Phi) %*% Y      # find maximum likelihood
}

#With the value of parameters w, and a given basis ??1,??2,. 
#we can estimate results for new points.

# compute estimates for points in x, given weigths W and respective basis
regression_curve <- function(xs, W, basisFn) {
  m <- sapply(basisFn,function(base) base(xs))  # values for each base
  apply(m, 1, function(row) sum(row*W))     # add them together, row by row
}

# function to draw regression line
draw_regression <- function(X, W, phi) {
  xs     <- seq(min(X),max(X),len=50)
  ys_hat <- regression_curve(xs, W, phi)
  points(xs, ys_hat, type="l", col="red")
}


#Let's define some basis functions
# some basis examples
one <- function(x) rep(1,length(x)) #phi-0 =1 (bias)
id  <- function(x) x
sq  <- function(x) x^2
x3  <- function(x) x^3
x4  <- function(x) x^4


#Standard Regression

# some data
X <- c(1,2,4,9,16)
Y <- c(1,9,144,293,921)

# basis for linear regression
phi <- c(one, id)
W <- compute_w(X, Y, phi)

plot(X,Y,pch=19)
abline(W, col="red")


#Polynomial regression. Notice that here, the basis are {??0,??1(x)=x,??2(x)=x^2}

# basis for quadratic regression
phi <- c(one, id, sq)
W <- compute_w(X, Y, phi)

plot(X,Y,pch=19)
draw_regression(X,W,phi)


# including a sine into the linear regression basis
phi <- c(one, id, function(x) sin(x))
W <- compute_w(X, Y, phi)

plot(X,Y,pch=19)
draw_regression(X,W,phi)
