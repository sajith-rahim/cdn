#PMRL Chapter 3

#----------------Linear Basis Models-----------------------
# Refer:
#    Sequentia Learning Gradient Descent
#----------------------------------------------------------
#----------------Linear Basis Models Sequential Learning---
update_w <- function(W, x_n, y_n, phi, eta=0.01) {
  phi_n <- matrix(sapply(phi, function(base) base(x_n)), ncol=1) # make it a vector so as to multiply
  W + eta * phi_n %*% (y_n - t(W)%*%phi_n)
}

compute_w_batch <-  function(X, Y, phi, convt=1e-3, eta=0.01) {
  W <- rnorm(length(phi),0,0.1)             # initialization to small random values
  for(i in 1:length(X)) {                   # batch update
    W <- update_w(W, X[i], Y[i], phi, eta)
  }
  W
}


X <- c(1,2,3,5,7)
Y <- c(3,5,6,12,21)
phi <- c(one, id) # basis for linear regression

W <- compute_w_batch(X, Y, phi,eta=0.015)
plot(X,Y,pch=19)
abline(W, col="red")
