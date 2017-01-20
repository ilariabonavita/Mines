# just testing

Cfunscript.path <- '/Users/ilaria_bonavita/My_project/tSNEcca/distr_cost_fun.R'
source(Cfunscript.path)

# compute derivative dP/dW1 (= dQ/dW2)
# every element is a matrix relative to ij distance



# compute derivative dC/dW1 (= dC/dW2)

deriv.3dmat <- function(Dist,W,allnum,distrvec){
  Bterm.3dmat <- array(rep(0,c(nrow(W)*ncol(W)*nrow(Dist))),c(nrow(W),ncol(W),nrow(Dist)))
  for(k in seq(nrow(Dist))){
    Bterm.3dmat[,,k] <- -2 * W%*%Dist[k,]%*%t(Dist[k,])*distrvec[k]*allnum[k] 
  }
  
  Bterm.mat <- apply(Bterm.3dmat,c(1,2), function(x) sum(x))
  
  deriv.3dmat <- array(rep(0,c(nrow(W)*ncol(W)*nrow(Dist))),c(nrow(W),ncol(W),nrow(Dist)))
  
  for(k in seq(nrow(Dist))){
    deriv.3dmat[,,k] <- Bterm.3dmat[,,k] - distrvec[k]*Bterm.mat
  }
  
  return(deriv.3dmat)
}

C.deriv <- function(Dist, W, allnum, distr1.vec, distr2.vec){
  
  
  deriv.3dmat <- function(Dist,W,allnum,distrvec){
  Bterm.3dmat <- array(rep(0,c(nrow(W)*ncol(W)*nrow(Dist))),c(nrow(W),ncol(W),nrow(Dist)))
  for(k in seq(nrow(Dist))){
    Bterm.3dmat[,,k] <- -2 * W%*%Dist[k,]%*%t(Dist[k,])*distrvec[k]*allnum[k] 
  }
  
  Bterm.mat <- apply(Bterm.3dmat,c(1,2), function(x) sum(x))
  
  deriv.3dmat <- array(rep(0,c(nrow(W)*ncol(W)*nrow(Dist))),c(nrow(W),ncol(W),nrow(Dist)))
  
  for(k in seq(nrow(Dist))){
    deriv.3dmat[,,k] <- Bterm.3dmat[,,k] - distrvec[k]*Bterm.mat
  }
  
  return(deriv.3dmat)
}
  distr1.deriv.3d <- deriv.3dmat(Dist,W,allnum,distr1.vec)
  d1overd2 <- distr1.vec/distr2.vec
  C <- matrix(0,dim(distr1.deriv.3d)[1],dim(distr1.deriv.3d)[2])
  for(k in seq(dim(distr1.deriv.3d)[3])){
  C <- C + 0.5*distr1.deriv.3d[,,k]*(log(d1overd2[k])+1-(d1overd2[k])^(-1))
  }
  return(C)
}

