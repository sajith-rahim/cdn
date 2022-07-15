#PMRL Chapter 8

#----------------Probabilistic Graphical Models------------
# Refer:
#    Probabilistic Graphical Models
#    Polynomial Curve Fitting 8.1.1
#----------------------------------------------------------
#--------Polynomial Curve Fitting -------------------------


library(BRugs)
library(coda)

#Sample Data
# sin(x) + Normal noise btw -pi to +pi
N  <- 50
xs <- seq(-pi,pi,len=N)
d  <- data.frame(x=xs,
                 y=sin(xs)+rnorm(N,0,0.1))
plot(d,pch=19)

data.list <- list(
  N  = N,
  alpha_1 = 10,
  sigma2 = 10,
  x = d$x,
  y = d$y
)


#Model

model = "
  model {
      for(i in 1:4) {
         w[i] ~ dnorm(0, alpha_1)       # prior for each w[i], w ~ N(0,1/alpha)
      }

      for(i in 1:N) {
          mu[i] <- w[1] + w[2]*x[i] + w[3]*pow(x[i],2) + w[4]*pow(x[i],3) 
          y[i] ~ dnorm(mu[i], sigma2)   # likelihood, y ~ N(mu, sigma^2)
      }
  }
"
#MCMC Params
samples=c("w")
chainLength=10000
burnin=0.10
n.chains=1
thin=1

#------BRugs --> BUGS : MCMC -------------------
  
writeLines(model, con="model.txt")  # write the model to a file
modelCheck( "model.txt" )           # Send the model to BUGS to check model syntax
modelData(bugsData(data.list))      
modelCompile(n.chains)              # BRugs command tells BUGS to compile the model


modelGenInits()                     # BRugs command tells BUGS to randomly initialize a chain


modelUpdate(chainLength*burnin)     # Burn-in period to be discarded
samplesSet(samples)                 # BRugs tells BUGS to keep a record of the sampled values
samplesSetThin(thin)                # Set thinning
modelUpdate(chainLength)            # BRugs command tells BUGS to randomly initialize a chain



# get posterior mean of p(w|data)
w_hat  <- samplesStats( "w" )$mean  
# to compute estimates for new values of y's based on a given w
compute_mu <- function(w,x) {
  w[1] + w[2]*x + w[3]*x^2 + w[4]*x^3
}
# for each x, estimate y given the posterior mean of p(w|data)
ys_hat <- sapply(xs, function(x) compute_mu(w_hat,x))

plot(d,pch=19)
points(xs,ys_hat,type="l", col="red",lwd=2)


#--We can produce a confidence interval, with the values of w computed by the mcmc,
#to get highest posterior density:
w_samples <- data.frame(w11=samplesSample("w[1]"))
for(i in 2:4)
  w_samples <- cbind(w_samples, samplesSample( paste0('w[',i,']') ))
names(w_samples) <- paste0("w",1:4)
head(w_samples)


plot(d,pch=19,type="n")
for(i in 1:20) {
  w <- w_samples[i,]
  y <- sapply(xs, function(x) compute_mu(w,x))
  points(xs,y,type="l", col="lightblue", lwd=1)
}
points(d,pch=19)

#Density region : 90% Confidence Interval

prob <- 0.9
hpd  <- matrix(rep(NA,N*2),ncol=2)
for(i in 1:N) {
  ys      <- apply(w_samples, 1, function(w) compute_mu(w,xs[i]))
  hpd[i,] <- HPDinterval(as.mcmc(ys), prob=prob)
}

plot(d,pch=19,type="n", xlab=paste0(100*prob,"% credible interval"), ylab="")
polygon(c(rev(xs), xs), c(rev(hpd[,1]), hpd[,2]), col = 'grey80', border = NA)
points(xs,ys_hat,type="l", col="red",lwd=2)
points(d,pch=19)
