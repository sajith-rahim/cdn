#Monte Carlo Markov Chain

# cf. Rizzo - Statistical Computing with R (2007)
metropolis.hastings <- function(f,  # the target distribution
                                g,  # the proposal distribution
                                rg, # a sample from the proposal distribution
                                x0, # initial value for chain, in R it is x[1]
                                chain.size=1e5,  # chain size
                                burn.perc=0.1) { # burn in percentage
  
  x <- c(x0, rep(NA,chain.size-1))  # initialize chain
  
  for(i in 2:chain.size)   {
    y <- rg(x[i-1])                 # generate Y from g(.|xt) using sampler rg
    alpha <- min(1, f(y)*g(x[i-1],y)/(f(x[i-1])*g(y,x[i-1])))
    x[i] <- x[i-1] + (y-x[i-1])*(runif(1)<alpha)  # update step
  }
  
  # remove initial part of the chain before output result
  x[(burn.perc*chain.size) : chain.size] 
}



a<-2.7; b<-6.3; size<-1e5

f  <- function(x)   dbeta(x,a,b)
rg <- function(x)   runif(1,0,1)
g  <- function(x,y) 1 #dunif(x,0,1) : for a unif(0,1) density = 1/(b-a)=1/1-0 

X <- metropolis.hastings(f,g,rg,x0=runif(1,0,1),chain.size=size)

par(mfrow=c(1,2),mar=c(2,2,1,1))
hist(X,breaks=50,col="blue",main="Metropolis-Hastings",freq=FALSE)
curve(dbeta(x,a,b),col="sienna",lwd=2,add=TRUE)
hist(rbeta(size,a,b),breaks=50,col="grey",main="Direct Sampling",freq=FALSE)
curve(dbeta(x,a,b),col="sienna",lwd=2,add=TRUE)


