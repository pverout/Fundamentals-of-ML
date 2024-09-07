set.seed(1)
X<-runif(60,-2,2) 
ee <- rnorm(60)
Y <- X^3 + ee #generating data
plot(X,Y)

activlogistic <- function(avec){ ##def of activation
  z <- (1+exp(-avec))^(-1)
  return(z)
}
deractivlogistic <- function(avec){ ##def of derivative of activation
  z <- exp(-avec)*(1+exp(-avec))^(-2)
  return(z)
}

nneteval <- function(xinvec,wmat1,wmat2,activfunc){ ##foward in the neural nework
  n.obs <- length(xinvec)
  X <- t(matrix(c(rep(1,n.obs),xinvec), nrow = n.obs, ncol = 2)) ##formating the matrix (adding a 1 for the  bias)
  a.1 <- wmat1%*%X #first layer a
  Z.1 <- activfunc(a.1) ##the hidden units in the first layer
  Z.1 <- rbind(rep(1,n.obs),Z.1) ##adding a 1 for the bias
  a.2 <- wmat2%*%Z.1 #second layer a
  return(a.2) #identify g(x) = x
}

nnetSSE <- function(yinc,xinvec,wmat1,wmat2,activfunc){
  y.hat <- as.vector(nneteval(xinvec,wmat1,wmat2,activlogistic)) # prediction made by neural net
  SSE <- sum((yinc - y.hat)^2) #squaring error
  return(SSE)
}


nnetGradient <- function(yinvec,xinvec,wmat1,wmat2,activfunc,deractivfunc){ ##finding gradient wrt to each element in w1, w2
  n.obs <- length(xinvec)
  X <- t(matrix(c(rep(1,n.obs),xinvec), nrow = n.obs, ncol = 2))
  a.1 <- wmat1%*%X ##needed for partials
  Z.1 <- activfunc(a.1)
  Z.1 <- rbind(rep(1,n.obs),Z.1) ##needed for partials
  
  sum1 <- matrix(0,nrow = dim(wmat1)[1], ncol = dim(wmat1)[2]) ##to compute the full gradient later
  sum2 <- matrix(0,nrow = dim(wmat2)[1], ncol = dim(wmat2)[2]) ##note the total gradient is just the sum of all observations gradient
  
  for(i in 1:length(xinvec)){ ##for loop over the observations
    del.O_del.g <- as.double(-2*(yinvec[i] - nneteval(xinvec[i],wmat1,wmat2,activfunc))) ##the first partial when taking the derivative of the O function
    del.g_del.aL <- 1 ##just to remind us that g'(x) = 1 since were using identity transformation
    
    #first layer
    del.aj1_del.wjk1 <- matrix(0,nrow = dim(wmat1)[1], ncol = dim(wmat1)[2]) ##partial when deriving the
    for(j in 1:dim(wmat1)[1]){ ##rows
      for(k in 1:dim(wmat1)[2]){ ##cols
        if(k == 1){
          del.aj1_del.wjk1[j,k] <- 1 ## in this case taking the derivative with respect to the bias is simply 1 
          next
        }
        else{
          del.aj1_del.wjk1[j,k] <- X[k,i] ## taking the derivative with respect to the weights yields in the first layer the X components
        }
      }
    }
    U2 <-1 ##recursive definition for U(2) initialized at 1
    
    #second layer
    del.aj2_del.wjk2 <- matrix(0,nrow = dim(wmat2)[1], ncol = dim(wmat2)[2])
    for(j in 1:dim(wmat2)[1]){ ##rows
      for(k in 1:dim(wmat2)[2]){ ##cols
        if(k == 1){ 
          del.aj2_del.wjk2[j,k] <- 1 ## in this case taking the derivative with respect to the bias is simply 1 
          next
        }
        else{
          del.aj2_del.wjk2[j,k] <- Z.1[k,i] ## taking the derivative with respect to the weights yields in the second layer the yields hidden units
        }
      }
    }
    U1 <- deractivfunc(a.1[,i])*wmat2[,2:6] ##recursive definition for U(s) keeping in mind U(2) is 1
    
    gradientW1 <- matrix(0,nrow = dim(wmat1)[1], ncol = dim(wmat1)[2])
    gradientW2 <- matrix(0,nrow = dim(wmat2)[1], ncol = dim(wmat2)[2])
    
    gradientW1[,1] <- U1*del.aj1_del.wjk1[,1] ##combining both partials to form (almost the whole gradient)
    gradientW1[,2] <- U1*del.aj1_del.wjk1[,2] ##combining both partials to form (almost the whole gradient)
    gradientW2 <- U2*del.aj2_del.wjk2
    
    gradientW11 <- del.O_del.g*gradientW1 ##multiplying finally by the derivative of the O function
    gradientW22 <- del.O_del.g*gradientW2
    
    sum1 <- sum1 + gradientW11 ##finding the sum with regards to each observation i to form the full gradient
    sum2 <- sum2 + gradientW22 ##gradient decomponsition
  }
  Gradient.list <- list(sum1,sum2) ##gradient output list
  return(Gradient.list)
}

set.seed(2)
w1 <- matrix(rnorm(10),nrow = 5,ncol = 2) #intializing our weights randomly
w2 <- matrix(rnorm(6),nrow = 1, ncol = 6) 
betas <- list(w1,w2)
nn <- 0.01 #learning rate
SSEvec <- c()
SSEvec[1] <- nnetSSE(Y,X,betas[[1]],betas[[2]],activlogistic)
for(J in 1:200){
  for(i in 1:length(X)){ ##iterating over data points
    delta.betas <- nnetGradient(Y[i],X[i],betas[[1]],betas[[2]],activlogistic,deractivlogistic) ##computing the gradient
    betas[[1]] <- betas[[1]] - nn*delta.betas[[1]] ##updating the betas based on their gradient
    betas[[2]] <- betas[[2]] - nn*delta.betas[[2]]
  }
  SSEvec[J+1] <- nnetSSE(Y,X,betas[[1]],betas[[2]],activlogistic) #recording the SSE after each epoch
}

plot(SSEvec, type = "l",
     xlab = "# Epochs", ylab = "SSE")
SSEvec

Y.hat <- nneteval(X,betas[[1]],betas[[2]],activlogistic)
plot(X, Y, col = "black",xlab = "x", ylab = "y", lwd = 2)
points(X, Y.hat, col = "blue", lwd = 2)
legend("topleft", legend = c("Actual", "Predicted"),
       col = c("black", "blue"), pch = 1, lwd = 2)

