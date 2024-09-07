set.seed(1)
n<-10000

#f functions
f1 <- function(x){
  return (4+3*sqrt(abs(x) + 2))
}
f2 <- function(x){
  return (sqrt(x + 12) + sin(x/2))
}
#initializing random portion
x1 <- runif(n,-7,7)
x2 <- runif(n,-12,12)
ee <- rnorm(n,0,1)
Y <- f1(x1) + f2(x2) + 0.5*ee
#ploting the 3D scatter and the f1 f2
dev.off()

## -- Function to be Estimated --
scatter3D(x1,x2,Y, theta =40, main = "y = f1(x) + f2(y)")
par(mfrow=c(2,1),  mar =c(3,3,1,1) + 1)
plot(x1,f1(x1))
plot(x2,f2(x2))

#creating an h function to sapply later
h <- function(x,ee){
  if(x<ee){
    return(0)
  }
  else{
    return ((x-ee)^3)
  }
}

## -- Builds the design matrix for cubic splines --
BuildSplines1DDesign <- function(xvec,nodes){
  DesignX <- matrix(0, nrow = length(xvec), ncol = length(nodes) + 4)
  #sum term without the h(x,c) function
  for(i in 0:3){
    DesignX[,i+1] <- t(xvec)^(i) 
  }
  #sum term with the h(x,c) function
  for(i in 4:(length(nodes) + 3)){
    DesignX[,i+1] <- t(sapply(xvec,h,ee=nodes[i-3])) #sapply makes evaluating the function for each node easy!
  }
  return (DesignX)
}
#intializing the nodes and the design matrices
nodesx1 <- seq(-5,5,by=2)
nodesx2 <- seq(-10,10,by=4)
Designx1 <- BuildSplines1DDesign(x1,nodesx1)
Designx2 <- BuildSplines1DDesign(x2,nodesx2)


#-- fhat outputs the y for splines --
fhat <- function(Xdesign,betas){
  y <- Xdesign%*%betas
  return(y)
}

##-- Backfitting procedure
betas1 <- rep(0,10)
betas2 <- rep(0,10)
for(i in 1:100){
  #first step
  U.l1 <- Y - fhat(Designx2,betas2)
  betas1 <- solve(t(Designx1)%*%Designx1)%*%t(Designx1)%*%U.l1
  
  #second step
  U.l2 <- Y - fhat(Designx1,betas1)
  betas2 <- solve(t(Designx2)%*%Designx2)%*%t(Designx2)%*%U.l2
}
betas1
betas2

xx1 <- seq(-7,7,by = 0.1)
xx2 <- seq(-12,12,by = 0.1)
#finding the estimated fhats by these two sequences above
Designxx1 <- BuildSplines1DDesign(xx1,nodesx1)
Designxx2 <- BuildSplines1DDesign(xx2,nodesx2)
f1hat <- fhat(Designxx1,betas1)
f2hat <- fhat(Designxx2,betas2)

#plots with correction constant
plot(xx1,f1(xx1), type = "l",lwd = 2, col = "black", xlab = "x1", ylab = "f(x1)", xlim=c(-7.5,7.5), ylim=c(7.5,17.5))
lines(xx1,f1hat, col = "blue", lwd = 3, lty = 5) 
title("Without Correction")
plot(xx2,f2(xx2), type = "l",lwd = 2, col = "black", xlab = "x2", ylab = "f(x2)", xlim=c(-12.5,12.5), ylim=c(-3,5))
lines(xx2,f2hat, col = "blue", lwd = 3, lty = 5)
legend(-12, 5, legend=c("Original", "Predicted"),
       col=c("black", "blue"), lwd = 2, lty=1:5, cex=0.5)


c1 <- mean(f1(xx1) - as.vector(Designxx1%*%betas1))
#plots without correction constant
plot(xx1,f1(xx1), type = "l",lwd = 2, col = "black", xlab = "x1", ylab = "f(x1)")
lines(xx1,f1hat + c1, col = "blue", lwd = 3, lty = 5) 
title("With Correction")
plot(xx2,f2(xx2), type = "l",lwd = 2, col = "black", xlab = "x2", ylab = "f(x2)")
lines(xx2,f2hat - c1, col = "blue", lwd = 3, lty = 5)
legend(-12, 5, legend=c("Original", "Predicted"),
       col=c("black", "blue"), lwd = 2, lty=1:5, cex=0.5)

