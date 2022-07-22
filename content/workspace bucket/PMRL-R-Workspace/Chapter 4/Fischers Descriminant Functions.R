#PMRL Chapter 4

#----------------Descriminant Functions-----------------------
# Refer:
#    Fischer's Linear Descriminant
#----------------------------------------------------------
#--------Fischer's Linear Descriminant---------------------


load("data.RData")

C_1 <- data[data[,3]==1,1:2]
C_2 <- data[data[,3]==2,1:2]

m_1 <- apply(C_1,2,mean)
m_2 <- apply(C_2,2,mean)

S_W <- matrix(rep(0,length(m_1)^2), ncol=2)

for(x_n in C_1)
  S_W <- S_W + (x_n - m_1) %*% t(x_n - m_1)
for(x_n in C_2)
  S_W <- S_W + (x_n - m_1) %*% t(x_n - m_1)

W <- solve(S_W) * (m_2-m_1)
W <- W / max(W)

#Project onto W
projection <- apply(data[,1:2],1,function(x) t(W)%*%x)


plot(0,0,type="n", xlim=c(1,6), ylim=c(-6,0))

points(projection[1,], projection[2,], col=data[,3], pch=3)

