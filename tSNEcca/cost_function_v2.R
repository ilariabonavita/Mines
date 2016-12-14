#### Define the cost function C


# Define pij 

p.fun.v2 <- function(W, dij, d, p.denom){
  p <- exp(-(norm(W%*%as.matrix(dij),'f')^2)/(2*sigma2))/p.denom
  return(p)
}

# Define qij
q.fun.v2 <- function(Z, dij, d, q.denom){
  q <- ((1 + norm(Z%*%as.matrix(dij),'f')^2)^(-1))/q.denom
  return(q)
}
  # define denominators
p.denom <- function(W, X.D){
  p.denom <- 0
  for(i in 1:nrow(X.D)){
    p.denom <- p.denom + exp(-(norm(W%*%as.matrix(X.D[i,]),'f')^2)/(2*sigma2))
  }
  return(p.denom)
}

q.denom <- function(Z, Y.D){
  q.denom <-0
  for (i in 1:nrow(Y.D)){
    q.denom <- q.denom + ((1+ norm(Z%*%Y.D[i,],'f')^2)^(-1))
  }
  return(q.denom)
}
  
# Define the cost function C
C.fun.v2 <- function(W,Z,X.D,Y.D, p.denom, q.denom){
  
  C <- 0
  for(i in 1:nrow(X.D)){
    C <- C+ p.fun.v2(W,X.D[i,],X.D,p.denom)*log(p.fun.v2(W,X.D[i,],X.D,p.denom)/q.fun.v2(Z,Y.D[i,],Y.D,q.denom))
  }
  return(C)
}


