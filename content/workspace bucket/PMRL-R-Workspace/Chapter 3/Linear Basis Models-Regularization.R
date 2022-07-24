#PMRL Chapter 3

#----------------Linear Basis Models-----------------------
# Refer:
    # Regularized Least Squares
    # Effect of Regularization L1/2,L1,L2...
#----------------------------------------------------------
#----------------Linear Basis Models Regularization--------

# compute weights for a given basis phi, using maximum likelihood
compute_w <- function(X, Y, phi) {
  Phi <- sapply(phi, function(basisFn) basisFn(X))  # make design matrix
  solve(t(Phi) %*% Phi) %*% t(Phi) %*% Y      # find maximum likelihood
}

# function to draw regression line
draw_regression <- function(X, W, phi) {
  xs     <- seq(min(X),max(X),len=50)
  ys_hat <- regression_curve(xs, W, phi)
  points(xs, ys_hat, type="l", col="red")
}

#----------------------------------------------------------


compute_w_reg <- function(X, Y, phi, lambda) {
  Phi <- sapply(phi, function(base) base(X))  # make design matrix
  solve(lambda * diag(length(phi)) + t(Phi) %*% Phi) %*% t(Phi) %*% Y  
  # w = [lambda*I + (phi)^T(phi)]^-1 * (phi)^T * Y
}

X <- c(1,2,3,5,7)
Y <- c(3,5,1,12,10)
phi <- c(one, id, sq, x3, x4)    # quartic regression

par(mfrow=c(1,2))

W <- compute_w(X, Y, phi)        # without regularization
plot(X,Y,pch=19,ylim=c(0,20), main="without regularization",sub="w =[(phi)^T(phi)]^-1 * (phi)^T * Y")
draw_regression(X,W,phi)

W <- compute_w_reg(X, Y, phi, lambda=11) # with regularization
plot(X,Y,pch=19,ylim=c(0,20), main="with L2 regularization",sub="w = [lambda*I + (phi)^T(phi)]^-1 * (phi)^T * Y")
draw_regression(X,W,phi)

