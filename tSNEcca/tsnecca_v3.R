library(entropy)


 N <- 10
 x <- runif(N) # latent signal
 r1 <- runif(N);  r2 <- runif(N)# random (helper) variables
 X <- cbind(tan(r1-x)+0.1*r1, r1+3*x-1/10*(sin(3*x)), r2*sin(x))
 Y <- cbind(x - 2*(1-exp(-r2))/(1+exp(-r2)), r2*x, tan(r2+x))

######################## check p function
X <- t(cbind(c(-1,2,1),c(2,5,1),c(-2,-1,4),c(2,3,1)))
X.pairs <- combn(seq(1:nrow(X)),2)
X.D <- X[X.pairs[1,],] - X[X.pairs[2,],]
W <- t(cbind( rep(1,3), 1:3))/4
pp <- exp(-(norm(W%*%(X[1,]-X[2,]))^2)/2)/(exp(-(norm(W%*%(X[1,]-X[2,]))^2)/2)+exp(-(norm(W%*%(X[1,]-X[3,]))^2)/2)+exp(-(norm(W%*%(X[2,]-X[3,]))^2)/2))
pp
pp.d <- (exp(-(norm(W%*%(X[1,]-X[2,]))^2)/2)+exp(-(norm(W%*%(X[1,]-X[3,]))^2)/2)+exp(-(norm(W%*%(X[2,]-X[3,]))^2)/2))
p.fun(W,(X[1,]-X[2,]))
p.deriv(W,(X[1,]-X[2,]))
#########################

#' @examples
NX <- 10
NY <- 10
q <- 3
# the dimension of the output transformed features
d <- 2

# generate two random views
X <- matrix(rnorm(NX*q),ncol=q)
Y <- matrix(rnorm(NY*q),ncol=q)


## remove the NAs in the matrices
if(anyNA(X)| anyNA(Y)){
  remove.na <- unlist(apply(cbind(X,Y), 1, anyNA))
  X <- X[!remove.na, , drop = F]
  Y <- Y[!remove.na, , drop = F]
  rm(remove.na) 
}

## Set default values to the unspecified parameters
# Random inizialisation

W <- matrix(rnorm(NX*q),nrow = d, ncol = ncol(X))
W <- matrix(0,nrow = d, ncol = ncol(X))
Z <- matrix(0,nrow = d, ncol = ncol(Y))
Z <- matrix(1,nrow = d, ncol = ncol(Y))

sigma2 <- 1
X.pairs <- combn(seq(1:nrow(X)),2)
Y.pairs <- combn(seq(1:nrow(Y)),2)

X.D <- X[X.pairs[1,],] - X[X.pairs[2,],]
Y.D <- Y[Y.pairs[1,],] - Y[Y.pairs[2,],]

# x adj
X.D <- X.D + max(X.D)

# Define the function to optimize
# 
p.fun <- function(W, X.d){ 
  
  p.denom <-sum( apply(X.pairs, 2, function(P){
    X.dkl <- X[P[1],]-X[P[2],]
    p.denom <- exp(-(norm(W%*%X.dkl)^2)/(2*sigma2))
    return(p.denom)
  }))
  p=exp(-(norm(W%*%X.d))^2/(2*sigma2))/p.denom
  return(p)
}





q.fun <- function(Z, Y.d){
  q.denom <- sum( apply(Y.pairs, 2, function(P){
    Y.dkl <- Y[P[1],] - Y[P[2],]
    q.denom <- (1 + norm(Z%*%Y.dkl)^2)^(-1)
    return(q.denom)
  }))
  q <- ((1 + norm(Z%*%Y.d)^2)^(-1))/q.denom
  return(q)
}

C.fun <- function(W,Z){
  C <- 0
  for(i in 1:ncol(X.pairs)){
    C <- C+ p.fun(W,X.D[,i])*log(p.fun(W,X.D[,i])/q.fun(Z,Y.D[,i]))
  }
  return(C)
}

p.deriv <- function(W,X.d){
  p1 <- W*0
  # p1 <-matrix(0,nrow = nrow(W), ncol = ncol(W))
  for(i in 1:ncol(X.pairs)){
    p1 <- - ( W %*% as.matrix(X.D[i,])) %*% t(as.matrix(X.D[i,]))*p.fun(W,X.D[i,]) + p1
  }
  p.deriv <- (-p.fun(W,X.d)/sigma2)*(W %*% as.matrix(X.d) %*% t(as.matrix(X.d)) + p1)
  return(p.deriv)
}

p.deriv.v2 <- function(W,X.d){
  p1 <- W*0
  # p1 <-matrix(0,nrow = nrow(W), ncol = ncol(W))
  for(i in 1:ncol(X.pairs)){
    p1 <- -2*W %*% as.matrix(X.D[i,]) %*% t(as.matrix(X.D[i,]))*p.fun(W,X.D[i,])*(p.fun(W,X.d)/sigma2) + p1
  }
  p.deriv <- (-p.fun(W,X.d)/sigma2)*(2*W %*% as.matrix(X.d) %*% t(as.matrix(X.d))) - p1
  return(p.deriv)
}

p.deriv.element <- function(W,X.d,t,z){
  p1 <- 0
  for (i in 1:ncol(X.pairs)){
    p1 <- p1 + (-2*W[t,]%*%X.D[i,]*(X.D[i,t])*p.fun(W,X.D[i,]))
  }
  p<- -(p.fun(W,X.d)/sigma2)*p1 - 2*W[t,]%*%X.d*(X.d[z])*p.fun(W,X.d)
  return(p)
}

q.deriv <- function(Z,Y.d){
  q1 <- Z*0
 # q1 <- matrix(0,nrow = nrow(Z), ncol = ncol(Z))
  for(i in 1:ncol(Y.pairs)){
    q1 <- q1 + (-q.fun(Z,Y.D[i,])*Z%*% as.matrix(Y.D[i,]) %*% t(as.matrix(Y.D[i,])))/(1+norm(Z%*%Y.D[i,])^2) 
  }
  q.deriv <- -2*q.fun(Z,Y.d)*( (Z%*%Y.d%*%t(Y.d))/(1+norm(Z%*%Y.d)^2) + q1)
  return(q.deriv)  
}


C.W.deriv <-  function(W,Z){
  C <- matrix(0,ncol = ncol(W), nrow = nrow(W))
  for (i in c(1:ncol(X.pairs))){
    C <- C+ p.deriv(W,X.D[,i])*(log(p.fun(W,X.D[,i])/q.fun(Z,Y.D[,i])) +1)
  }
  return(sum(C))
}

C.Z.deriv <- function(W,Z){
  C <- matrix(0,ncol = ncol(Z), nrow = nrow(Z))
  for (i in c(1:ncol(Y.pairs))){
    C <- C+ (q.deriv(Z,Y.D[,i])/q.fun(Z,Y.D[,i]))*p.fun(W,X.D[,i])
  }
  return(sum(C))
}

#=========================================
# gradient checking
set.seed(124)

W.test <- matrix(rnorm(NX*q),nrow = d, ncol = ncol(X))
W.test <- t(apply(W.test,1, function(x){ x/sum(x)}))
Z.test <- matrix(rnorm(NY*q),nrow = d, ncol = ncol(Y))
Z.test <- t(apply(Z.test,1, function(x){ x/sum(x)}))
Z.test <- matrix(1,nrow = d, ncol = ncol(Y))
Z.test[2,2] <- 2.0


W.test <- matrix(1,nrow = d, ncol = ncol(X))
W.test  <- matrix(sample.int(4, d*ncol(X), replace = T), nrow = d, ncol = ncol(X))
eps <-0.0001

t <- 1
z <- 1
W.plus <- W.test 
W.plus[t,z] <- W.test[t,z ] + eps
W.minus <- W.test 
W.minus[t,z] <- W.test[t,z] - eps

Z.plus <- Z.test 
Z.plus[t,z] <- Z.test[t,z] + eps
Z.minus <- Z.test 
Z.minus[t,z] <- Z.test[t,z] - eps


X.test <- X.D[1,]
Y.test <- Y.D[1,]


W.test <- W
# p.deriv
p.deriv.approx <-(p.fun(W.plus,X.test) - p.fun(W.minus, X.test))/(2*eps) 
p.deriv(W.test,X.test)
p.deriv.approx

par(new=TRUE)
plot(k,p.deriv(W.test,X.test)[t,z])
par(new=TRUE)
plot( k, p.deriv.approx, type="l", col="green")


# q.deriv
q.deriv.approx <-(q.fun(Z.plus,Y.test) - q.fun(Z.minus, Y.test))/(2*eps) 
q.deriv(Z.test,Y.test)
q.deriv.approx

