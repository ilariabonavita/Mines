##### Compute the derivatives 
source('/Users/ilaria_bonavita/Mines/tSNEcca/cost_function_v2.R')

p1.deriv <- function(W,d,p.denom){
p1 <- W*0
for(i in 1:nrow(d)){
  p1 <- -W %*% as.matrix(d[i,]) %*% t(as.matrix(d[i,]))*p.fun.v2(W,d[i,],d, p.denom) + p1
}
return(p1)
}

q1.deriv <- function(Z, d, q.denom){
  q1 <- Z*0
  for (i in 1:nrow(d)){
    q1 =  q1 + (-q.fun.v2(Z,d[i,],d,q.denom)*Z %*% d[i,] %*% t(d[i,]))/(1+ norm(Z%*%d[i,],'f')^2)
  }
return(q1)
}

# pij derivative
p.deriv.v2 <- function(W,dij, d, p.denom, p1){
  p.deriv <- (p.fun.v2(W,dij,d, p.denom)/sigma2)*(-W %*% as.matrix(dij) %*% t(as.matrix(dij)) - p1)
  return(p.deriv)
}

# qij derivative
q.deriv.v2 <- function(Z,dij,d,q.denom, q1){
  q.deriv <- (-2*q.fun.v2(Z,dij,d, q.denom))*((Z%*%dij%*%t(dij))/(1+ norm(Z%*%dij,'f')^2) + q1)
  return(q.deriv)  
}

# C derivative with respect to W
C.W.deriv.v2 <-  function(W,Z,dx,dy,p.denom,q.denom,p1){
  C <- W*0
  for (i in c(1:nrow(dx))){
    C <- C+ p.deriv.v2(W,dx[i,],dx,p.denom,p1)*(log(p.fun.v2(W,dx[i,],dx,p.denom)/q.fun.v2(Z,dy[i,],dy,q.denom)) +1)
  }
  return(C)
}

# C derivative with respect to Z
C.Z.deriv.v2 <- function(W,Z,dx,dy,p.denom,q.denom,q1){
  C <- W*0
  for (i in c(1:nrow(dx))){
    C <- C+ (-q.deriv.v2(Z,dy[i,],dy,q.denom,q1)/q.fun.v2(Z,dy[i,],dy,q.denom))*p.fun.v2(W,dx[i,],dx,p.denom)
  }
  return(C)
}

#=========================================
# # gradient checking
# set.seed(1234)
# 
# W.test <- matrix(rnorm(NX*q),nrow = d, ncol = ncol(X))*10
# W.test <- t(apply(W.test,1, function(x){ x/sum(x)}))
# Z.test <- matrix(rnorm(NY*q),nrow = d, ncol = ncol(Y))
# Z.test <- t(apply(Z.test,1, function(x){ x/sum(x)}))
# 
# W.test <- matrix(2,nrow = d, ncol = ncol(X))
# W.test  <- matrix(sample.int(4, d*ncol(X), replace = T), nrow = d, ncol = ncol(X))
# eps <-0.0001
# 
# t <- 1
# z <- 1
# W.test <- W
# W.plus <- W.test 
# W.plus[t,z] <- W.test[t,z] + eps
# W.minus <- W.test 
# W.minus[t,z] <- W.test[t,z] - eps
# 
# Z.test <- Z
# Z.plus <- Z.test 
# Z.plus[t,z] <- Z.test[t,z] + eps
# Z.minus <- Z.test 
# Z.minus[t,z] <- Z.test[t,z] - eps
# 
# 
# X.test <- X.D[12,]
# Y.test <- Y.D[12,]
# 
# # p.deriv
# p.deriv.approx <-(p.fun(W.plus,X.test,X.D) - p.fun(W.minus, X.test,X.D))/(2*eps) 
# p.deriv(W.test,X.test,X.D)[t,z]
# p.deriv.approx
# 
# par(new=TRUE)
# plot(k,p.deriv(W.test,X.test)[t,z])
# par(new=TRUE)
# plot( k, p.deriv.approx, type="l", col="green")
# 
# 
# # q.deriv
# q.deriv.approx <-(q.fun(Z.plus,Y.test, Y.D) - q.fun(Z.minus, Y.test, Y.D))/(2*eps) 
# q.deriv(Z.test,Y.test, Y.D)
# q.deriv.approx
# 
# # C.W deriv
# C.W.deriv.approx <- (C.fun(W.plus, Z, X.D, Y.D)- C.fun(W.minus, Z, X.D, Y.D))/(2*eps) 
# C.W.deriv.approx
# C.W.deriv(W.test,Z,X.D,Y.D)
# 
# # C.Z deriv
# C.Z.deriv.approx <- (C.fun(W, Z.plus, X.D, Y.D)- C.fun(W, Z.minus, X.D, Y.D))/(2*eps) 
# C.Z.deriv.approx
# C.Z.deriv(W,Z.test,X.D,Y.D)


