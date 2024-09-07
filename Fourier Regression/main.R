setwd("C:\\Users\\pvero\\OneDrive\\Desktop\\Concordia\\Winter 2024\\FunOfML\\Fourier Regression")
load("Loaddata.RData")
str(Loadseries)
Weekseries
str(Weekseries)
length(Weekseries)

##-- Initializing the Design Matrix --##
P <-5
Ct.sin.j <- function(t,j){sin(3*pi/2 + 2*pi*j*t/(365.25/7))}
Ct.cos.j <- function(t,j){cos(3*pi/2 + 2*pi*j*t/(365.25/7))}
#Initializing Design Matrix
Designmatrix <- matrix(0,nrow = length(Weekseries), ncol=2*P+1)
#Loading the first Column with 1s
Designmatrix[,1] <- 1
#Loading the design Matrix (approach: [1, sin.js,cos.js])
for(t in 0:length(Weekseries)-1){
  for(j in 1:P){
    Designmatrix[t+1,2*j-1 + 1] <- Ct.sin.j(t,j)
  }
  for(j in 1:P){
    Designmatrix[t+1,2*j + 1] <- Ct.cos.j(t,j)
  }
}
Designmatrix


##-- Creating the OLSE Funtion --##
OLSSE <- function(yTrain,XdesignTrain,XdesignValid,yValid){
  betas <- solve(t(XdesignTrain)%*%XdesignTrain)%*%t(XdesignTrain)%*%yTrain
  y.hat <- XdesignValid%*%betas
  SSE <- t(y.hat-yValid)%*%(y.hat-yValid)
  return(SSE)
}
