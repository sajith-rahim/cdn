#PMRL Chapter 4

#----------------Probabilistic Generative Models-----------
# Refer:
#    Probabilistic Generative Models
#    MLE
#----------------------------------------------------------
#--------Probabilistic Generative Models ------------------
load("data.RData")


sigmoid <- function(x) { 1/(1+exp(-x)) }
curve(sigmoid, -5, 5, col="red", lwd=2)

get_posteriors <- function(x, priors, mus, Sigma) {
  K <- length(priors)
  S_inv <- solve(Sigma)
  
  a_k <- rep(NA,K)
  for (k in 1:K) {
    w_k    <- S_inv %*% mus[k,]
    w_k0   <- -0.5 * t(mus[k,]) %*% S_inv %*% mus[k,] + log(priors[k])
    a_k[k] <- t(w_k) %*% x + w_k0
  }
  sum_posteriors <- sum(sapply(1:K, function(k) exp(a_k[k])))
  sapply(1:K, function(k) exp(a_k[k]))/sum_posteriors
}


mle_mus <- function(X, Y) {
  N_k <- table(Y)
  K   <- length(N_k)
  
  t(sapply(1:K, function(k) apply(X[Y==k,],2,sum)/N_k[k]))
}

mle_Sigma <- function(X, Y, mus) {
  N_k <- table(Y)
  K   <- length(N_k)
  N   <- nrow(X)
  D   <- ncol(X)
  Sigma <- matrix(rep(0,D^2), D)
  
  for(k in 1:K) {                  # for a given class C_k
    S_k <- matrix(rep(0,D^2), D) 
    X_k <- X[Y==k,]                # select all x's from class C_k
    for(i in 1:nrow(X_k)) {        # for each x in C_k
      S_k <- S_k + (X_k[i,] - mus[k,]) %*% t(X_k[i,] - mus[k,])
    }
    Sigma <- Sigma + S_k
  }
  
  Sigma/N
}

X <- data[,1:2]
Y <- data[,3]

# get maximum likelihood estimates
mus   <- mle_mus(X,Y)
Sigma <- mle_Sigma(X,Y,mus)

# define equal priors for C_1 and C_2
priors <- c(0.5,0.5)

# get a classifier function for this dataset
classifier <- function(x) { get_posteriors(x, priors, mus, Sigma) }

# create some new data and classify them - here 3 data points
X_new <- matrix(c(4.0,5.0,  
                  6.0,3.0,  
                  5.3,4.5), ncol=2, byrow=T)

round(t(apply(X_new, 1, classifier)), 2)


