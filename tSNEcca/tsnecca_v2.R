# tSNE-cca method
# description Compute tSNE-correlation analysis between two data matrices
# 
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

##### Set default values to the unspecified parameters

# Random inizialisation
W <- matrix(rnorm(NX*q),nrow = d, ncol = ncol(X))
Z <- matrix(rnorm(NY*q),nrow = d, ncol = ncol(Y))

# Normilize matrices
W <- t(apply(W,1, function(x){ x/sum(x)}))
Z <- t(apply(Z,1, function(x){ x/sum(x)}))

# Assign value to sigma^2 (to do: tune parameter)
sigma2 <- 1

#### Define a matrix of the distance between Xi and Xj (Yi and Yj) for all the possible combination of elements in X (Y)
X.pairs <- combn(seq(1:nrow(X)),2)
Y.pairs <- combn(seq(1:nrow(Y)),2)

# For simplicity, we use only the distances from now on
X.D <- X[X.pairs[1,],] - X[X.pairs[2,],]
Y.D <- Y[Y.pairs[1,],] - Y[Y.pairs[2,],]

#### Define the cost function C

# Define pij 
p.fun <- function(W, dij, d){
  p1 <- 0
  for(i in 1:nrow(d)){
    p1 <- p1 + exp(-(norm(W%*%as.matrix(d[i,]),'f')^2)/(2*sigma2))
  }
  p <- exp(-(norm(W%*%as.matrix(dij),'f')^2)/(2*sigma2))/p1
  return(p)
}
 
# Define qij
q.fun <- function(Z, dij, d){
  q.denom <-0
  for (i in 1:nrow(d)){
    q.denom <- q.denom + ((1+ norm(Z%*%d[i,],'f')^2)^(-1))
  }
  q <- ((1 + norm(Z%*%as.matrix(dij),'f')^2)^(-1))/q.denom
  return(q)
}

# Define the cost function C
C.fun <- function(W,Z,X.D,Y.D){
  C <- 0
  for(i in 1:nrow(X.D)){
    C <- C+ p.fun(W,X.D[i,],X.D)*log(p.fun(W,X.D[i,],X.D)/q.fun(Z,Y.D[i,],Y.D))
  }
  return(C)
}

##### Compute the derivatives 

# pij derivative
p.deriv <- function(W,dij, d){
  p1 <- W*0
  for(i in 1:nrow(d)){
    p1 <- -W %*% as.matrix(d[i,]) %*% t(as.matrix(d[i,]))*p.fun(W,d[i,],d) + p1
    }
  p.deriv <- (p.fun(W,dij,d)/sigma2)*(-W %*% as.matrix(dij) %*% t(as.matrix(dij)) - p1)
  return(p.deriv)
}

# qij derivative
q.deriv <- function(Z,dij,d){
  q1 <- Z*0
  for (i in 1:nrow(d)){
    q1 =  q1 + (-q.fun(Z,d[i,],d)*Z %*% d[i,] %*% t(d[i,]))/(1+ norm(Z%*%d[i,],'f')^2)
  }
  q.deriv <- (-2*q.fun(Z,dij,d))*((Z%*%dij%*%t(dij))/(1+ norm(Z%*%dij,'f')^2) + q1)
return(q.deriv)  
}

# C derivative with respect to W
C.W.deriv <-  function(W,Z,dx,dy){
  C <- W*0
  for (i in c(1:nrow(dx))){
    C <- C+ p.deriv(W,dx[i,],dx)*(log(p.fun(W,dx[i,],dx)/q.fun(Z,dy[i,],dy)) +1)
  }
  return(C)
}

# C derivative with respect to Z
C.Z.deriv <- function(W,Z,dx,dy){
  C <- W*0
  for (i in c(1:nrow(dx))){
    C <- C+ (-q.deriv(Z,dy[i,],dy)/q.fun(Z,dy[i,],dy))*p.fun(W,dx[i,],dx)
  }
  return(C)
}

#=========================================
# gradient checking
set.seed(1234+k)

W.test <- matrix(rnorm(NX*q),nrow = d, ncol = ncol(X))*10
W.test <- t(apply(W.test,1, function(x){ x/sum(x)}))
Z.test <- matrix(rnorm(NY*q),nrow = d, ncol = ncol(Y))
Z.test <- t(apply(Z.test,1, function(x){ x/sum(x)}))

W.test <- matrix(2,nrow = d, ncol = ncol(X))
W.test  <- matrix(sample.int(4, d*ncol(X), replace = T), nrow = d, ncol = ncol(X))
eps <-0.0001

t <- 1
z <- 1
W.test <- W
W.plus <- W.test 
W.plus[t,z] <- W.test[t,z] + eps
W.minus <- W.test 
W.minus[t,z] <- W.test[t,z] - eps

Z.test <- Z
Z.plus <- Z.test 
Z.plus[t,z] <- Z.test[t,z] + eps
Z.minus <- Z.test 
Z.minus[t,z] <- Z.test[t,z] - eps


X.test <- X.D[12,]
Y.test <- Y.D[1,]

# p.deriv
p.deriv.approx <-(p.fun(W.plus,X.test) - p.fun(W.minus, X.test))/(2*eps) 
p.deriv(W.test,X.test)[t,z]
p.deriv.approx

par(new=TRUE)
plot(k,p.deriv(W.test,X.test)[t,z])
par(new=TRUE)
plot( k, p.deriv.approx, type="l", col="green")


# q.deriv
q.deriv.approx <-(q.fun(Z.plus,Y.test, Y.D) - q.fun(Z.minus, Y.test, Y.D))/(2*eps) 
q.deriv(Z.test,Y.test, Y.D)
q.deriv.approx

# C.W deriv
C.W.deriv.approx <- (C.fun(W.plus, Z, X.D, Y.D)- C.fun(W.minus, Z, X.D, Y.D))/(2*eps) 
C.W.deriv.approx
C.W.deriv(W.test,Z,X.D,Y.D)

# C.Z deriv
C.Z.deriv.approx <- (C.fun(W, Z.plus, X.D, Y.D)- C.fun(W, Z.minus, X.D, Y.D))/(2*eps) 
C.Z.deriv.approx
C.Z.deriv(W,Z.test,X.D,Y.D)



#################### check ok