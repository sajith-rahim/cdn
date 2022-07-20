#---------------Rejection Sampling-------------------


# generate n samples from f using rejection sampling with g
accept.reject <- function(f, c, g, rg, n) { 
  n.accepts     <- 0
  result.sample <- rep(NA, n)
  
  while (n.accepts < n) {
    y <- rg(1)               # Generate a y from bounding dist. g
    u <- runif(1,0,1)        # Generate random uniform
    if (u < f(y)/(c*g(y))) { # Accept or RejectA
      n.accepts <- n.accepts+1
      result.sample[n.accepts] = y
    }
  }
  
  result.sample
}



#----------------Sampling p(x) = Beta(2,2)--------------------
# p(x) = Beta(2,2)
# q(x) = U
# C = 2
f  <- function(x) 6*x*(1-x)     # pdf of Beta(2,2), maximum density is 1.5
g  <- function(x) x/x           # g(x) = 1 but in vectorized version
rg <- function(n) runif(n,0,1)  # uniform, in this case
c  <- 2                         # c=2 since f(x) <= 2 g(x)

vals <- accept.reject(f, c, g, rg, 10000) 

# Plotting
hist(vals, breaks=30, freq=FALSE, main="Sample vs true Density")
xs <- seq(0, 1, len=100)
lines(xs, dbeta(xs,2,2), col="red", lwd=2)


#--------------Plot Acceptance-Rejection------------------------

# xs <- seq(0, 1, len=100)
# plot(xs, dbeta(xs,2,2), ylim=c(0,c*1.3), type="l", col="red", lwd=2, ylab="densities")
# lines(xs, c*g(xs), type="l", col="blue", lwd=2)
# legend("topleft",c("f(x)","c*g(x)"), col=c("red","blue"), lwd=2) 
# 
# draw.segment <- function(begin.segment, end.segment) {
#   segments(c(begin.segment,end.segment,end.segment,begin.segment), c(0,0,c*1.025,c*1.025), 
#            c(end.segment,end.segment,begin.segment,begin.segment), c(0,c*1.025,c*1.025,0))
#   n.pts <- 1000
#   us <- runif(n.pts, 0, 1)
#   ys <- begin.segment + rg(n.pts)*(end.segment-begin.segment)
#   accepted <- us < f(ys)/(c*g(ys))
#   points(ys, c*us, col=ifelse(accepted,"green","red"), pch=19)  
# }
# 
# draw.segment(0.00, 1.0) 

#c=10 ---------------Poor acceptance ratio-------------------------------------

# c <- 10
# 
# xs <- seq(0, 1, len=100)
# plot(xs, dbeta(xs,2,2), ylim=c(0,c*1.25), type="l", col="red", lwd=2, ylab="densities")
# lines(xs, c*g(xs), type="l", col="blue", lwd=2)
# legend("topleft",c("f(x)","c*g(x)"), col=c("red","blue"), lwd=2) 
# 
# draw.segment(0.10, 0.20)
# draw.segment(0.45, 0.55)
# draw.segment(0.90, 1.00)
