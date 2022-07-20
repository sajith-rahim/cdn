# MCMC - Hastings Algorithm
# 
# Target Distribution = weighted sum of two normal distributions. 
# 
# This sort of distribution is fairly straightforward to sample from, 
# but let's draw samples with MCMC. 
# 
# This is a contrived example, but distributions like this are not totally impossible, 
# and might arise when sampling things from a mixture 
# (such as human heights, which are bimodal due to sexual dimorphism).
# 
# Fairly arbitrarily, here are 
# parameters and the definition of the target density.


p <- 0.4
mu <- c(-1, 2)
sd <- c(.5, 2)
#MIXTURE
f <- function(x)
  p     * dnorm(x, mu[1], sd[1]) +
  (1-p) * dnorm(x, mu[2], sd[2])

curve(f(x), col="red", -4, 8, n=301, las=1)


#proposal :  normal distribution centred on the current point
#with a standard deviation of 4

q <- function(x) rnorm(1, x, 4)

#Metropolis Algorithm
step <- function(x, f, q) {
  ## Pick new point
  xp <- q(x)
  ## Acceptance probability:
  alpha <- min(1, f(xp) / f(x))
  ## Accept new point with probability alpha:
  if (runif(1) < alpha)
    x <- xp
  ## Returning the point:
  x
}


run <- function(x, f, q, nsteps) {
  res <- matrix(NA, nsteps, length(x))
  for (i in seq_len(nsteps))
    res[i,] <- x <- step(x, f, q)
  drop(res)
}

res <- run(-10, f, q, 1000)
layout(matrix(c(1, 2), 1, 2), widths=c(4, 1))
par(mar=c(4.1, .5, .5, .5), oma=c(0, 4.1, 0, 0))
plot(res, type="s", xpd=NA, ylab="Parameter", xlab="Sample", las=1)
usr <- par("usr")
xx <- seq(usr[3], usr[4], length=301)
plot(f(xx), xx, type="l", yaxs="i", axes=FALSE, xlab="")


hist(res, 50, freq=FALSE, main="", ylim=c(0, .4), las=1,
     xlab="x", ylab="Probability density")
z <- integrate(f, -Inf, Inf)$value
curve(f(x) / z, add=TRUE, col="red", n=200)



set.seed(1)
res.long <- run(-10, f, q, 50000)
hist(res.long, 100, freq=FALSE, main="", ylim=c(0, .4), las=1,
     xlab="x", ylab="Probability density", col="grey")
z <- integrate(f, -Inf, Inf)$value
curve(f(x) / z, add=TRUE, col="red", n=200)


#100, 1,000, 10,000 and 100,000 steps

n <- 10^(2:5)
samples <- lapply(n, function(n) run(-10, f, q, n))
xlim <- range(sapply(samples, range))
br <- seq(xlim[1], xlim[2], length=100)

hh <- lapply(samples, function(x) hist(x, br, plot=FALSE))
ylim <- c(0, max(f(xx)))

par(mfrow=c(2,2), mar=rep(.5, 4), oma=c(4, 4, 0, 0))
for (h in hh) {
  plot(h, main="", freq=FALSE, yaxt="n",
       ylim=range(h$density, ylim))
  curve(f(x), add=TRUE, col="red", n=300)
}






#--------------------2D Gaussian---------------------
# fn makes a multivariate normal density given a 
# vector of means (centre of the distribution) 
# and variance-covariance matrix.



make.mvn <- function(mean, vcv) {
  logdet <- as.numeric(determinant(vcv, TRUE)$modulus)
  tmp <- length(mean) * log(2 * pi) + logdet
  vcv.i <- solve(vcv)
  
  function(x) {
    dx <- x - mean
    exp(-(tmp + rowSums((dx %*% vcv.i) * dx))/2)
  }
}


#Mixture Unweigted

mu1 <- c(-1, 1)
mu2 <- c(2, -2)
vcv1 <- matrix(c(1, .25, .25, 1.5), 2, 2)
vcv2 <- matrix(c(2, -.5, -.5, 2), 2, 2)
f1 <- make.mvn(mu1, vcv1)
f2 <- make.mvn(mu2, vcv2)
f <- function(x)
  f1(x) + f2(x)

x <- seq(-5, 6, length=71)
y <- seq(-7, 6, length=61)
xy <- expand.grid(x=x, y=y)
z <- matrix(apply(as.matrix(xy), 1, f), length(x), length(y))

image(x, y, z, las=1)
contour(x, y, z, add=TRUE)



#Assume that we don't actually know how to sample from a mvn 
# (it's not actually hard, but this is simpler), 
# let's make a proposal distribution that is uniform in two dimensions, 
# sampling from the square with width 'd' on each side.


q <- function(x, d=8)
  x + runif(length(x), -d/2, d/2)

x0 <- c(-4, -4)
set.seed(1)
samples <- run(x0, f, q, 1000)

image(x, y, z, xlim=range(x, samples[,1]), ylim=range(x, samples[,2]))
contour(x, y, z, add=TRUE)
lines(samples[,1], samples[,2], col="#00000088")


#100000 samples
samples <- run(x0, f, q, 100000)


smoothScatter(samples)
#contour(x, y, z, add=TRUE)