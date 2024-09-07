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


## -- K-Fold cross Val --##
SSEmat <- matrix(data = 0, nrow = 5, ncol = 6)
dates <- c("2008-01-01", "2009-01-01", "2010-01-01",
           "2011-01-01", "2012-01-01","2013-01-01")
cutoffs <- c()
cutoffs[1] <- 1
for(i in 2:7){
  cutoffs[i] <- max(which(Weekseries < dates[i-1]))
}
cutoffs 
##attempted to use cuttoffs to find week on which to prefrom k- fold cross validation


for(i in 1:5){ ## hard coded k fold cross validation (the i's in this for loop iterate over potetial P values)
  SubDesignmatrix <- Designmatrix[,1:(2*i+1)]
  SSEmat[i,1]<-OLSSE(Loadseries[54:291], SubDesignmatrix[54:291,], SubDesignmatrix[1:53,],Loadseries[1:53] ) 
  SSEmat[i,2]<-OLSSE(Loadseries[c(1:53,106:291)], SubDesignmatrix[c(1:53,106:291),],SubDesignmatrix[54:105,] ,Loadseries[54:105])
  SSEmat[i,3]<-OLSSE(Loadseries[c(1:105,158:291)], SubDesignmatrix[c(1:105,158:291),],SubDesignmatrix[106:157,] ,Loadseries[106:157])          
  SSEmat[i,4]<-OLSSE(Loadseries[c(1:157,210:291)], SubDesignmatrix[c(1:157,210:291),],SubDesignmatrix[158:209,] ,Loadseries[158:209])
  SSEmat[i,5]<-OLSSE(Loadseries[c(1:209,262:291)], SubDesignmatrix[c(1:209,262:291),],SubDesignmatrix[210:261,] ,Loadseries[210:261])
  SSEmat[i,6]<-OLSSE(Loadseries[c(1:261)], SubDesignmatrix[c(1:261),],SubDesignmatrix[262:291,] ,Loadseries[262:291])
}
SSEmat
SSEtot <- rowSums(SSEmat)
SSEtot

plot(SSEtot, xlab = "Model Complexity",main = "SSE over Fourier Space")
which.min(SSEtot)


##-- Plotting preds --##
P.op <- which.min(SSEtot)
Designmatrix
Designmatrix.op <- Designmatrix[,1:(2*P.op+1)]
betas <- solve(t(Designmatrix.op)%*%Designmatrix.op)%*%t(Designmatrix.op)%*%Loadseries
OptPreds <- Designmatrix.op%*%betas
dev.off()
##plot the prediction and emprical values
plot(Loadseries,col = "black", type = "l", lwd = 2, xlab = "Enumerated Weeks", ylab = "Energy Consumption")
lines(OptPreds, col = "blue", lwd = 3, lty = 5)
legend(1, 7600000, legend=c("Original", "Predicted"),
       col=c("black", "blue"), lwd = 2, lty=1:5, cex=0.8)
