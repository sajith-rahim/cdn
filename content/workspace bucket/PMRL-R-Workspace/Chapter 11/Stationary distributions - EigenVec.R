#--------------Stationary distributions---------------


# transition matrix for 3 states
#P[i,j] gives the probability of moving from state i to state j
P <- rbind(c(.5,  .25, .25),
           c(.2,  .1,  .7),
           c(.25, .25, .5))

# #rows sum up to one i.e the chain must go somewhere,
# even if that place is the same place. 
rowSums(P) 
#not 
colSums(P)


# This function takes a state vector x 
# (where x[i] is the probability of being in state i) and 
# iterates it by multiplying by the transition matrix P, 
# advancing the system for n steps.

iterate.P <- function(x, P, n) {
  res <- matrix(NA, n+1, length(x))
  res[1,] <- x
  for (i in seq_len(n))
    res[i+1,] <- x <- x %*% P
  res
}

n <- 10
#start from state 1 i.e [1,0,0] 100%prob of in state 1
y1 <- iterate.P(c(1, 0, 0), P, n)
# Similarly for other states
y2 <- iterate.P(c(0, 1, 0), P, n)
y3 <- iterate.P(c(0, 0, 1), P, n)


matplot(0:n, y1, type="l", lty=1, xlab="Step", ylab="y", las=1)
matlines(0:n, y2, lty=2)
matlines(0:n, y3, lty=3)



# We can use R's eigen function to extract the leading eigenvector for the system
# (the t() here transposes the matrix so that we get the left eigenvector).
v <- eigen(t(P), FALSE)$vectors[,1]
v <- v/sum(v) # normalise eigenvector

#Lets plot to see how cloase we are
matplot(0:n, y1, type="l", lty=1, xlab="Step", ylab="y", las=1)
matlines(0:n, y2, lty=2)
matlines(0:n, y3, lty=3)
points(rep(10, 3), v, col=1:3)


#----------------------------------------------

# The proceedure above iterated the overall probabilities of different states;
# not the actual transitions through the system. So, let's iterate the system,
# rather than the probability vector. The function run here takes a state
# (this time, just an integer indicating which of the states $1, 2, 3$ the system is in),
# the same transition matrix as above, and a number of steps to run.
# Each step, it looks at the possible places that it could transition to and
# chooses 1 (this uses R's sample function).


run <- function(i, P, n) {
  res <- integer(n)
  for (t in seq_len(n))
    res[[t]] <- i <- sample(nrow(P), 1, pr=P[i,])
  res
}

samples <- run(1, P, 1000)
plot(samples, type="s", xlab="Step", ylab="State", las=1)

cummean <- function(x)
  cumsum(x) / seq_along(x)

plot(cummean(samples == 1), type="l", ylim=c(0, 1),
     xlab="Step", ylab="y", las=1)
lines(cummean(samples == 2), col=2)
lines(cummean(samples == 3), col=3)
