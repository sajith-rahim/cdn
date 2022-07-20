#---------------Importance Sampling-------------------


i.sampling <- function(f, g, h, rh, n=1e4) {
  ys <- rh(n)
  mean(g(ys)*f(ys)/h(ys))
}


#P(X>4.5)



xs <- seq(0.1,10,by=0.05)
plot(xs,dexp(xs),col="blue", type="l")   # the exponential pdf
lines(xs,exp(-(xs-4.5)),col="red",lwd=2) # the truncated pdf
abline(v=4.5,lty=2)