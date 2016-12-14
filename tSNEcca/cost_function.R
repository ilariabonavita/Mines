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
