#### Define the cost function C

#--- W, Z just for testing
W <- matrix(rnorm(N*q),nrow = d, ncol = ncol(XTr))
W <- t(apply(W,1, function(x){ x/sum(x)}))
Z <- matrix(rnorm(N*q),nrow = d, ncol = ncol(YTr))
Z <- t(apply(Z,1, function(x){ x/sum(x)}))

X.D <- XTr.D 
Y.D <- YTr.D

n <- 10
#----------------------------

# Define pij and qij

### things that can be computed only once

WXD <- rep(0,nrow(X.D))
WXD <- apply(X.D, 1, function(x){
  wdij <- norm(W%*%as.matrix(x),'f')^2
  return(wdij)
})

ZYD <- rep(0,nrow(Y.D))
ZYD <- apply(Y.D, 1, function(x){
  zdij <- norm(Z%*%as.matrix(x),'f')^2
  return(zdij)
})

p.numer <- apply(as.array(WXD),1, function(i){
  exp(-i/(2*sigma2))})

q.numer <- apply(as.array(ZYD),1,function(i){
  (1+i)^(-1)
  })
q.denom <- sum(q.numer)

### Define pij
pij.fun <- function(WXD, i, j, p.numer, XTr.pairs,n){
  p.num <- p.numer[which((XTr.pairs[1,]==i & XTr.pairs[2,]==j) | (XTr.pairs[1,]==j & XTr.pairs[2,]==i))]
  pji.denom <- sum(as.array(p.numer)[which(XTr.pairs[1,]==i | XTr.pairs[2,]==i)])
  pij.denom <- sum(as.array(p.numer)[which(XTr.pairs[1,]==j | XTr.pairs[2,]==j)])
  pjcondi <- p.num/pji.denom
  picondj <- p.num/pij.denom
  pij <- (pjcondi + picondj)/(2*n)
  return(pij)
  }

### Define qij
qij.fun <- function(ZYD,i,j,q.numer,q.denom){
 qij <- q.numer[which((YTr.pairs[1,]==i & YTr.pairs[2,]==j) | (YTr.pairs[1,]==j & YTr.pairs[2,]==i))]/q.denom
 return(qij)
 }

### Define the cost function C
C.fun <- function(){
  C <- 0
  for(i in 1:ncol(XTr.pairs)){
    C <- C+ pij.fun(WXD,XTr.pairs[1,i],XTr.pairs[2,i],p.numer,XTr.pairs,n) *log(pij.fun(WXD,XTr.pairs[1,i],XTr.pairs[2,i],p.numer,XTr.pairs,n)/qij.fun(ZYD,YTr.pairs[1,i],YTr.pairs[2,i],q.numer,q.denom))
  }
  return(C)
}
