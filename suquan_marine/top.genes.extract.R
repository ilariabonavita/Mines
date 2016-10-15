# feature selection

top.genes.extract <- function(xtrain,ytrain, ngenes){

plist <- lapply(1:ncol(xtrain),function(i){
  fit <-lm(xtrain[,i]~ytrain)
  f <- summary(fit)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail = F)
  attributes(p) <- NULL
  return(p)}
  )

pvec <- unlist(plist)

idx.tokeep <- order(pmat)[1:ngenes]
return(idx.tokeep)
}
