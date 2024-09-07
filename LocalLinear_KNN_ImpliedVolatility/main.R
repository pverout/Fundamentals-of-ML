setwd("C:\\Users\\pvero\\OneDrive\\Desktop\\Concordia\\Winter 2024\\FunOfML\\LocalLinear_KNN_ImpliedVolatility")
load("IV.RData")
library(plot3D)
library(matlib)

##-- Gaussian Kernel Function --##
dmvn <- function(xmat,muvec,sigmamat){
  p <- length(muvec)
  pdfvec <- c()
  for (i in 1:nrow(xmat)){
    pdfvec[i] <- exp(-0.5*(xmat[i,]-muvec)%*%inv(sigmamat)%*%t(xmat[i,]-muvec))/(sqrt(((2*pi)^p)*det(sigmamat)))
  }
  return(pdfvec)
}

##-- Local Linear Regression Betas --
gridMaturity <- seq(from = 0, to = 1, by = 0.05)
gridStrike <- seq(from = 850, to = 1150, by = 20)
##OLSSE.weighted using the formula derived from the exercises
OLSSE.weighted <- function(X,Y,W){
  betas <- inv(t(X)%*%W%*%X)%*%t(X)%*%W%*%Y
  return(betas)
}

##storing the covariance matrices for later
Cov.mat <- list(
  matrix(c(0.3,0,0,100), nrow = 2, ncol = 2),
  matrix(c(0.1,0,0,25), nrow = 2, ncol = 2),
  matrix(c(0.9,0,0,200), nrow = 2, ncol =2))
fhatmatlist <- list()
Xdesign <- matrix(c(rep(1,192),Maturity, Strike), nrow = 192, ncol = 3)
muvec <- matrix(c(mean(Maturity), mean(Strike)),ncol = 2, nrow = 1)
Cov.mat

##-- Local Linear Regression on Cov matrices -- 
for(i in 1:3){
  sigmamat <- Cov.mat[[i]]
  tempmat <- matrix(0,nrow = 21, ncol = 16)
  for(j in 1:21){ ##iterates through the set of maturities
    for(k in 1:16){ ##iterates through the set of strikes
      muvec <- matrix(c(gridMaturity[j],gridStrike[k]), nrow = 1, ncol = 2) #computing the mu vector
      W <- diag(dmvn(Xdesign[,-1],muvec,sigmamat)) # computing weights
      betas <- OLSSE.weighted(Xdesign,ImplicitVol,W) # preforming local linear regression
      tempmat[j,k] <- betas[1] + betas[2]*gridMaturity[j] + betas[3]*gridStrike[k] ##making the prediction
    }
  }
  fhatmatlist[[i]] <- tempmat ##storing the predictions in the right matrix
}##Runs long (shouldn't be over 10 seconds)...Complexity O(n^4)





##-- Plotting Results of Local Linear Regression --
length(gridMaturity)
##here we need to make matrices for the x and y coordiates for the surf function
MatrixMaturity <- matrix(gridMaturity, nrow = 21, ncol = 16)
MatrixMaturity
MatrixStrike <- t(matrix(gridStrike, nrow =16, ncol = 21))
MatrixStrike

par(mfrow=c(2,2), oma = c(1,0,0,1) + 0.1, mar =c(2,0,1,1) + 0.3)
surf3D(MatrixStrike,MatrixMaturity,fhatmatlist[[1]], ticktype = "detailed", cex.axis = 0.5, cex.lab = 0.5, bty = "b2", ylab = "Time to Maturity", xlab = "Strike", zlab = "Implied Volatility", theta = 40)
title("Sigma 1")
surf3D(MatrixStrike,MatrixMaturity,fhatmatlist[[2]], ticktype = "detailed", cex.axis = 0.5, cex.lab = 0.5, bty = "b2", ylab = "Time to Maturity", xlab = "Strike", zlab = "Implied Volatility", theta = 40)
title("Sigma 2")
surf3D(MatrixStrike,MatrixMaturity,fhatmatlist[[3]], ticktype = "detailed", cex.axis = 0.5, cex.lab = 0.5, bty = "b2", ylab = "Time to Maturity", xlab = "Strike", zlab = "Implied Volatility", theta = 40)
title("Sigma 3")
scatter3D(Strike,Maturity,ImplicitVol, ticktype = "detailed", cex.axis = 0.5, cex.lab = 0.5, bty = "b2", ylab = "Time to Maturity", xlab = "Strike", zlab = "Implied Volatility", theta = 40)
title("Historical Data")



##-- Defining a Nearest Neighbor Function --

##distance function
distance <-  function(a,b){return(sum((a-b)^2)^(1/2))}
##returns the indices of K nearest neighbours
KNN.set <- function(x,obs,k) {
  dist <- apply(x,1, distance,obs)
  NN.ind <- which(dist %in% sort(dist)[1:k])
  return(NN.ind)
}
##preforms predictions based on the set of K nearest neighbours
PredictNearestNeighbors <- function(yval,xval,xpred,NNnumber){
  NNpred <- c()
  for(i in 1:dim(xpred)[1]){
    NN.ind <- KNN.set(xval,xpred[i,],NNnumber)
    NNpred[i] <- mean(yval[NN.ind])
  }
  return(NNpred)
}

##-- Preforming KNN for different Ks -- 
xpred.matrix <- matrix(c(rep(gridMaturity, each = 16),rep(gridStrike, 21)), nrow = 21*16, ncol = 2)
xpred.matrix
Xdesign <- matrix(c(rep(1,192),Maturity, Strike), nrow = 192, ncol = 3)
Xdesign[,3] <-  Xdesign[,3]/400 ## dividing strike by 400 to normalize
xpred.matrix[,2] <- xpred.matrix[,2]/400

fhatmatlistKNN <- list()
Ks <- c(5,15,30)
for(i in 1:3){
  KNNpred.vec <- PredictNearestNeighbors(ImplicitVol,Xdesign[,-1],xpred.matrix,Ks[i]) ##preforms predictions
  fhatmatlistKNN[[i]] <- matrix(KNNpred.vec,nrow  = 21, ncol = 16, byrow=TRUE)
  
}


##-- Plotting KNN Results
par(mfrow=c(2,2), oma = c(1,0,0,1) + 0.1, mar =c(2,0,1,1) + 0.3)
surf3D(MatrixStrike,MatrixMaturity,fhatmatlistKNN[[1]], ticktype = "detailed", cex.axis = 0.5, cex.lab = 0.5, bty = "b2", ylab = "Time to Maturity", xlab = "Strike", zlab = "Implied Volatility", theta = 40)
title("K = 5")
surf3D(MatrixStrike,MatrixMaturity,fhatmatlistKNN[[2]], ticktype = "detailed", cex.axis = 0.5, cex.lab = 0.5, bty = "b2", ylab = "Time to Maturity", xlab = "Strike", zlab = "Implied Volatility", theta = 40)
title("K = 15")
surf3D(MatrixStrike,MatrixMaturity,fhatmatlistKNN[[3]], ticktype = "detailed", cex.axis = 0.5, cex.lab = 0.5, bty = "b2", ylab = "Time to Maturity", xlab = "Strike", zlab = "Implied Volatility", theta = 40)
title("K = 20")
scatter3D(Strike,Maturity,ImplicitVol, ticktype = "detailed", cex.axis = 0.5, cex.lab = 0.5, bty = "b2", ylab = "Time to Maturity", xlab = "Strike", zlab = "Implied Volatility", theta = 40)
title("Historical Data")

