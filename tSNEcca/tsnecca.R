
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
 W <- matrix(,nrow = d, ncol = ncol(X))
 W <- matrix(1,nrow = d, ncol = ncol(X))
 Z <- matrix(,nrow = d, ncol = ncol(Y))
 Z <- matrix(1,nrow = d, ncol = ncol(Y))
 
 sigma2 <- 1
 X.pairs <- combn(seq(1:nrow(X)),2)
 Y.pairs <- combn(seq(1:nrow(Y)),2)
 
 X.D <- apply(X.pairs, 2, function(P){
       X.dkl <- X[P[1],]-X[P[2],]
       return(X.dkl)})
 Y.D <- apply(Y.pairs, 2, function(P){
   Y.dkl <- Y[P[1],]-Y[P[2],]
   return(Y.dkl)})
 
 # Define the function to optimize
# 
p.fun <- function(W, X.d){ 
    p.denom <-sum( apply(X.pairs, 2, function(P){
    X.dkl <- X[P[1],]-X[P[2],]
    p.denom <- exp(-norm(t(W)*X.dkl)^2/(2*sigma2))
    return(p.denom)
    }))   
    p=exp(-norm(t(W)*X.d)^2/(2*sigma2))/p.denom
 return(p)
 }
# 
#    
# P.mat.fun <- function(p.fun){
#   P <- apply(X.pairs, 2, function(P){
#   pp <- p(W, X[P[1],]-X[P[2],],p.denom)
#   return(pp)
#   })
#   return(P)}
# 
# 
# q.fun <- function(Z, Y.d){
#   q.denom <- sum( apply(Y.pairs, 2, function(P){
#   Y.dkl <- Y[P[1],]-Y[P[2],]
#   q.denom <- exp(-norm(t(Z)*Y.dkl)^2)
#   return(q.denom)
# }))  
# q=exp(-norm(t(Z)*Y.d)^2)/q.denom
# return(q)
# }
# 
# Q.mat.fun <- function(q.fun){
#   Q <- apply(Y.pairs, 2, function(P){
#     qq <- q(Z, Y[P[1],]-Y[P[2],])
#     return(qq)
#   })
#   return(Q)}


P <- function(W){
  P.num <- apply(X.pairs, 2, function(x){
    X.dij <- X[x[1],] - X[x[2],]
    pp <- exp((-norm(t(W)*X.dij)^2)/(2*sigma2))
    return(pp)
  })
  P.den <- sum(P.num)
  
  return(P.num/P.den)
}

Q <- function(Z){
  Q.den <- apply(Y.pairs, 2, function(y){
    Y.dij <- Y[y[1],] - Y[y[2],]
    qq <- (1+norm(t(Z)*Y.dij)^2)
    return(qq)
  })
  Q.num <- sum(Q.den)
  
  return(Q.num/Q.den)
}


C.fun <- function(W,Z){
  
 sum(P(W)*(log(P(W)/Q(Z))))
  
}

P.deriv <- function(W){
  
  p1 <-matrix(0,nrow = nrow(W), ncol = ncol(W))
  for(i in 1:ncol(X.pairs)){
    p1 <- p1 + (-W %*% as.matrix(X.D[,i]) %*% t(as.matrix(X.D[,i])))*p.fun(W,D[,i])
  }
  
  p.deriv <- matrix(0,nrow = nrow(W), ncol = ncol(W))
  for(i in 1:ncol(X.pairs)){
    p.deriv <- (p.fun(W,X.D[,i])/sigma2)*(-W %*% as.matrix(X.D[,i]) %*% t(as.matrix(X.D[,i])) - p1)
  }
  
  return(p.deriv)
}
  
  
C.W.deriv <- function(W,Z){
 sum(P.deriv(W)%*%(log(P(W)/Q(Z)))+1)
}
 

P.deriv(W)
 