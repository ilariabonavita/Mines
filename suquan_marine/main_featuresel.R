library(ROCR)
library(glmnet)
library(suquan)
library(apg)
library(parallel)  #need to load!
set.seed(1234)

mc.cores <- 3
data.path <- '/Users/ilaria_bonavita/Mines/mines_yunlong/data/arloth.txt'
setwd('/Users/ilaria_bonavita/Mines/suquan/')
source('./R/predict.suquan.R')
source('./top.genes.extract.R')

exp.path <- paste0(getwd(),'/feature_selection/test1/')

#--- load data
dataset <- read.table(data.path,h=T)
x <- dataset[,-1]
y <- dataset[,1]
n <- nrow(x)
p <- ncol(x)

#--- Make train/test splits
nrepeats <- 10
nfolds <- 3
folds <- list()
for (i in seq(nrepeats)) {
  folds <- c(folds,split(sample(seq(n)), rep(1:nfolds, length = n)))
}
nexp <- length(folds)

#--- Set lambdas parameters
#lambda <- 1.5^{-10:20}/(p*n*(nfolds-1)/nfolds)
#lambda <- 1.5^{-seq(18:40)}/(p*n*(nfolds-1)/nfolds)
#nlambda <- length(lambda)

lambda <- 1.307067e-07
niter <- 2

#--- set number of genes 
genes <- seq(10,300,20)
ngenes <- length(genes)


#--- save parameters for reproducibility 
param <- list(nrepeats=nrepeats, nfolds = nfolds, lambda= lambda, genes=genes)
saveRDS(param,paste0(exp.path,'parameters.RDS'))
#---


# Main loop (done in parallel on several cores)
perf <- mclapply(seq(nexp), function(iexp){
  auc <- matrix(nrow=nlambda, ncol=niter)
  acc <- matrix(nrow=nlambda, ncol=niter)
  itest <- folds[[iexp]]
  itrain <- seq(n)[-itest]
  for (iiter in seq(niter)) {
    # iter is how many times we optimize the quantile distribution
    
    for (iigenes in seq(ngenes)) {
      # lambda is the regularization parameter
      tokeep <- top.genes.extract(x[itrain,],y[itrain],genes[iigenes])
      # Train SUQUAN on the training set
      m <- suquan(x[itrain,tokeep],y[itrain], family="binomial",opts=list(alpha=0), lambda=lambda, maxiter=iiter)
      
      # Predict on the test set
      ypred <- predict.suquan(m,x[itest,tokeep])
      
      # Compute AUC on the test set
      pred <- prediction(ypred, y[itest])
      auc[ilambda, iiter] <- performance(pred, "auc")@y.values[[1]]
      
      # Compute accuracy on the test set
      acc.all <- performance(pred, "acc")@y.values[[1]]
      acc.max <- acc.all[which.max(acc.all)]
      acc[ilambda,iiter] <- acc.max
    }
  }
  
  return(list(auc=auc,acc=acc))
},mc.cores=mc.cores)

acc.list <- list()
auc.list <- list()
for (i in seq(nexp)){
  perf.exp <- perf[[i]]
  acc.list[[i]] <- list(perf.exp$acc)
  auc.list[[i]] <- list(perf.exp$auc)
}

saveRDS(acc.list,paste0(exp.path,'acc_list.RDS'))
saveRDS(auc.list,paste0(exp.path,'auc_list.RDS'))

mauc <- apply(array(unlist(auc.list), dim = c(nrow(auc.list[[1]][[1]]), ncol(auc.list[[1]][[1]]), length(auc.list))), c(1,2), mean)
pdf(paste0(exp.path,"auc.pdf"), width=5, height=5)
matplot(log(lambda), mauc, type="l", lty=1, lwd=2, main="MDD", ylab="AUC",xlab="log(lambda)")
legend("topright", legend=c("Gaussian", "SUQUAN"), col=seq(niter), lty=1, lwd=2)
dev.off()

macc <- apply(array(unlist(acc.list), dim = c(nrow(acc.list[[1]][[1]]), ncol(acc.list[[1]][[1]]), length(acc.list))), c(1,2), mean)
pdf(paste0(exp.path,"acc.pdf"), width=5, height=5)
matplot(log(lambda), macc, type="l", lty=1, lwd=2, main="MDD", ylab="acc",xlab="log(lambda)")
legend("topright", legend=c("Gaussian", "SUQUAN"), col=seq(niter), lty=1, lwd=2)
dev.off()

