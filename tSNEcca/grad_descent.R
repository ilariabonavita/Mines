# libraries
library(sgd)

# Parameters
gamma <- 0.01
maxiter <- 2
precision <- 0.1

# Initialize W and Z (and normalize)
W <- matrix(rnorm(NX*q),nrow = d, ncol = ncol(X))
W <- t(apply(W,1, function(x){ x/sum(x)}))
Z <- matrix(rnorm(NY*q),nrow = d, ncol = ncol(Y))
Z <- t(apply(Z,1, function(x){ x/sum(x)}))

Wold <- W
Zold <- Z
C.fun(Wold,Zold)
for(i in c(1: maxiter)){
  W <- Wold - gamma*C.W.deriv(Wold,Zold)
  Z <- Zold - gamma*C.Z.deriv(Wold,Zold)
  print(W)
  print(Z)
  Wold <- W
  Zold <- Z
}

C.fun(W,Z)
C.fun(Wold,Zold)

W1 <- apply(W,1, function(x){ x/sum(x)})
