#### Define the cost function C

# Define pij and qij

### things that can be computed only once

### Define Student t distribution

ij.distr <- function(Dist, W, id.pair, allnumer.vec, denom){
  
  distr <- allnumer.vec[id.pair]/denom
  
  return(distr)
}

allnumerator.vec <- function(Dist, W){
  numerator <- rep(0,nrow(Dist))
  numerator <- apply(Dist, 1, function(x){
    numij <- (1 + norm(W%*%as.matrix(x),'f')^2)^(-1)
    return(numij)
  })
}


distr.vec <- function(Dist,W,allnum){
  denom <- sum(allnum)
  d <- rep(0,nrow(Dist))
  for (i in seq(nrow(Dist))){
    d[i] <- ij.distr(Dist,W,i,allnum,denom)
  }
  return(d)
}

## Define the cost function
C.fun <- function(P.distr,Q.distr){

KL.fun <- function(P.distr, Q.distr){
  KL <- 0
  for(i in 1:length(P.distr)){
    KL <- KL+ P.distr[i]*log(P.distr[i]/Q.distr[i])
  }
  return(KL)
}### first compute KL(P||Q)
KLpq <- KL.fun(P.distr, Q.distr)
KLqp <- KL.fun(Q.distr, P.distr)
  C<- 0.5*(KLpq + KLqp)
  return(C)
}

# #================================================================================
# #--- example for testing only
# 
# ## Set some parameters for data generation
# N <- 10 # overall number of examples
# MaxAngle <- 4*pi
# MinRadius <-0.3
# MaxRadius <- 8
# sx <- 0.5
# sy <- 0.5
# set.seed(1234)
# ## Generate data for views 1,2
# t <- seq(0, MaxAngle, length.out = N)
# r <- seq(MinRadius, MaxRadius, length.out = N) + 2*runif(N)
# #### generate X, the noise can be added!
# X <- cbind(r*cos(t+0*rnorm(N)*0.05), r*sin(t+0*rnorm(N)*0.05),rnorm(N*2))
# #### generate Y, the noise can be added!
# Y <- cbind(t+0*rnorm(N)*1, 2*rnorm(N),rnorm(N*2))
# ## Generate two initial random projcetion matrices
# W <- matrix(rnorm(N*q),nrow = d, ncol = ncol(XTr))
# W <- t(apply(W,1, function(x){ x/sum(x)}))
# Z <- matrix(rnorm(N*q),nrow = d, ncol = ncol(YTr))
# Z <- t(apply(Z,1, function(x){ x/sum(x)}))
# 
# ## Define distance matrices for the two views
# X.pairs <- combn(seq(1:nrow(X)),2)
# Y.pairs <- combn(seq(1:nrow(Y)),2)
# X.Dist <-  X[X.pairs[1,],] - X[X.pairs[2,],]
# Y.Dist <- Y[Y.pairs[1,],] - Y[Y.pairs[2,],]
# 
# ## Compute distribution of pairwise distances in the projected space
# X.distr.vec <- distr.vec(X.Dist,W)
# Y.distr.vec <- distr.vec(Y.Dist,Z)
# 
# ## Compute Cost function
# # compute the two KL matrices
# KLpq <- KL.fun(X.distr.vec, Y.distr.vec)
# KLqp <- KL.fun(Y.distr.vec, X.distr.vec)
# 
# # compute the cost function
# C <- C.fun(KLpq,KLqp)
# 
# 
