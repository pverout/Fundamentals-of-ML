setwd("C:\\Users\\pvero\\OneDrive\\Desktop\\Concordia\\Winter 2024\\FunOfML\\YieldCurve PCA")
load("YieldCurvedata.RData")


length(Maturitiesvec)
YCmatrix<-data.matrix(YCdf[,-1]) ##setting up the design matrix
Allmeans <- colMeans(YCmatrix)
YCmatrixdm <- scale(YCmatrix, scale = FALSE) ##demeaning the data without rescalling

pca <- princomp(YCmatrixdm) ##applying PCA
pca$sdev[1:6] ##standard devition of each component
cum.exp.var <- cumsum(pca$sdev[1:6]^2)/sum(pca$sdev^2)  #cumulative variance
cum.exp.var

a1 <- pca$loadings[,1] #extracting the loadings
a2 <- pca$loadings[,2]
a3 <- pca$loadings[,3]

dev.off()
plot(1:120, a1, type = "l",lwd = 2, col = "orange", xlab = "time to maturity", ylab = "yield", ylim=c(-0.12,0.2))
lines(1:120,a2, col = "blue", lwd = 2) 
lines(1:120, a3, col = "red", lwd = 2, lty = 2)
legend(0.1,0.2,legend=c("1st Loading", "2nd Loading", "3rd Loading"),
       col=c("orange", "blue", "red"), lwd = 2, lty=1:5, cex=0.5)

data.matrix(pca$loadings)
X <- t(t(data.matrix(pca$loadings[,1:120]))%*%t(YCmatrixdm))
score.matrix <- pca$scores
score.matrix
index1 <- which(YCdf[,1] == "1991-01-02") ##finiding the index at which dates
index2 <- which(YCdf[,1] == "2000-05-31")
index3 <- which(YCdf[,1] == "2003-03-26")
index4 <- which(YCdf[,1] == "2006-06-27")


pca.X <- function(index,p){ ##function that computes x based on only p components 
  x <- numeric(120)
  for(i in 1:p){
    x <- x + pca$loadings[,i]*score.matrix[index ,i] ##computes the x based on the loading vecs and their respective scores
  }
  return(x+Allmeans)
}
##pca for 3 components
pca.X(index1,3)
pca.X(index2,3)
pca.X(index3,3)
pca.X(index4,3)

#check using all 120 components
sum(pca.X(index1,120)-YCmatrix[index1,]) #practically zero
sum(pca.X(index2,120)-YCmatrix[index2,])
sum(pca.X(index3,120)-YCmatrix[index3,])
sum(pca.X(index4,120)-YCmatrix[index4,])

##using 3 principal components
par(mfrow=c(2,2))
plot(1:120, pca.X(index1,3), type = "l",lwd = 2, col = "black", xlab = "time to maturity", ylab = "yield", ylim = c(0.097,0.106), main = "1991-01-02")
lines(1:120,YCmatrix[index1,], col = "blue", lwd = 2, lty = 5) 
plot(1:120, pca.X(index2,3), type = "l",lwd = 2, col = "black", xlab = "time to maturity", ylab = "yield", ylim = c(0.056,0.062), main = "2000-05-31")
lines(1:120,YCmatrix[index2,], col = "blue", lwd = 2, lty = 5)
plot(1:120, pca.X(index3,3), type = "l",lwd = 2, col = "black", xlab = "time to maturity", ylab = "yield", main = "2003-03-26")
lines(1:120,YCmatrix[index3,], col = "blue", lwd = 2, lty = 5) 
plot(1:120, pca.X(index4,3), type = "l",lwd = 2, col = "black", xlab = "time to maturity", ylab = "yield", main = "2006-06-27")
lines(1:120,YCmatrix[index4,], col = "blue", lwd = 2, lty = 5)
text(0.5,0.5,"Second title",cex=2,font=2)
mtext("3-PCA", side = 3, line = -1.5, outer = TRUE)

##using 5 principal components
par(mfrow=c(2,2))
plot(1:120, pca.X(index1,5), type = "l",lwd = 2, col = "black", xlab = "time to maturity", ylab = "yield", ylim = c(0.097,0.106), main = "1991-01-02")
lines(1:120,YCmatrix[index1,], col = "blue", lwd = 2, lty = 5) 
plot(1:120, pca.X(index2,5), type = "l",lwd = 2, col = "black", xlab = "time to maturity", ylab = "yield", ylim = c(0.056,0.062), main = "2000-05-31")
lines(1:120,YCmatrix[index2,], col = "blue", lwd = 2, lty = 5)
plot(1:120, pca.X(index3,5), type = "l",lwd = 2, col = "black", xlab = "time to maturity", ylab = "yield", main = "2003-03-26")
lines(1:120,YCmatrix[index3,], col = "blue", lwd = 2, lty = 5) 
plot(1:120, pca.X(index4,5), type = "l",lwd = 2, col = "black", xlab = "time to maturity", ylab = "yield", main = "2006-06-27")
lines(1:120,YCmatrix[index4,], col = "blue", lwd = 2, lty = 5)
text(0.5,0.5,"Second title",cex=2,font=2)
mtext("5-PCA", side = 3, line = -1.5, outer = TRUE)





