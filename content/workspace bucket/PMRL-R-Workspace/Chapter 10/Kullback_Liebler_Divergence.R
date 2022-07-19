#----------------Variational Inference------------
# Refer:
#    10
#    
#----------------------------------------------------------
#--------KL Divergence-------------------------------------

# KL for continuous functions
KL <- function(q, p, lower, upper, ...) {
  f <- function(x) q(x, ...) * log(q(x, ...)/p(x))
  integrate(f, lower, upper)$value
}

# an eg:
p <- function(x) dgamma(x, shape=3.0, scale=0.25)
q1 <- function(x) dlnorm(x, 0, 1)
q2 <- function(x) dlnorm(x, 0, .45)

curve(p, 0, 6, lwd=2, ylab="")
curve(q1,      0, 6, lwd=2, col="red",   add=T)
curve(q2,      0, 6, lwd=2, col="green", add=T)  




KL(q1,p, lower=1e-3, upper=100)
## [1] 1.709245
KL(q2,p, lower=1e-3, upper=100)
## [1] 0.3400462
