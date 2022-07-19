#----------------Variational Inference------------
# Refer:
#    10
#    
#----------------------------------------------------------
#-------------Variational Inference------------------------
#--Example: Approximate a gamma using log-normals----------

# Kullback Liebler Divergence
KL <- function(q, p, lower, upper, ...) {
  f <- function(x) q(x, ...) * log(q(x, ...)/p(x))
  integrate(f, lower, upper)$value #Sum-over for discrete
}


#-----------------ALGORITHM-----------------------------
#1. Select a family of distributions Q
#2. Find best approximation of posterior unnormalized p say q(z) as
#           arg-min qEQ KL(q(z) || p(z) )
#

variational_lnorm <- function(p_unnormalized, lower, upper) {
  q <- dlnorm # Let q be a log-normal
  
  J <- function(params) {
    KL(q, p_unnormalized, lower=lower, upper=upper, meanlog=params[1], sdlog=params[2])
  }
  
  optim(par=c(0, 1), fn=J)$par
}





p_unnormalized <- function(x) dgamma(x, shape=3.0, scale=0.25)

approximation_params <- variational_lnorm(p_unnormalized, lower=1e-3, upper=100)


# approximated distribution:
q <- function(x) dlnorm(x, approximation_params[1], approximation_params[2])



KL(q,p_unnormalized,1e-3,10) # KL 
## [1] 0.02765858



curve(p_unnormalized, 0, 6, lwd=2, ylab="", ylim=c(0,1.25))
curve(q,       0, 6, lwd=2, col="red",   add=T)
