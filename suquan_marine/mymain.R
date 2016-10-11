library(ROCR)
library(glmnet)
library(suquan)
set.seed(1234)
mc.cores <- 3

# Download some data
source("https://bioconductor.org/biocLite.R")
biocLite("breastCancerTRANSBIG")
library(breastCancerTRANSBIG)
data(transbig)
x <- t(exprs(transbig))

# Label is 10-year metastasis
dataset <- read.table('/Users/ilaria_bonavita/Mines/mines_yunlong/data/arloth.txt')

x <- x[,-1]
y <- y[,1]
n <- nrow(x)
p <- ncol(x)

# Make train/test splits
nrepeats <- 10
nfolds <- 3
folds <- list()
for (i in seq(nrepeats)) {
  folds <- c(folds,split(sample(seq(n)), rep(1:nfolds, length = n)))
}
nexp <- length(folds)

lambda <- 1.5^{-10:20}/(p*n*(nfolds-1)/nfolds)
nlambda <- length(lambda)
niter <- 2

# Main loop (done in parallel on several cores)
auc <- mclapply(seq(nexp), function(iexp){
  auc <- matrix(nrow=nlambda, ncol=niter)
  itest <- folds[[iexp]]
  itrain <- seq(n)[-itest]
  for (iiter in seq(niter)) {
    # iter is how many times we optimize the quantile distribution

    for (ilambda in seq(nlambda)) {
      # lambda is the regularization parameter

      # Train SUQUAN on the training set
      m <- suquan(x[itrain,],y[itrain], family="binomial",opts=list(alpha=0), lambda=lambda[ilambda], maxiter=iiter)

      # Predict on the test set
      ypred <- predict.suquan(m,x[itest,])

      # Compute AUC on the test set
      pred <- prediction(ypred, y[itest])
      auc[ilambda, iiter] <- performance(pred, "auc")@y.values[[1]]
    }
  }
  return(auc)
},
mc.cores=mc.cores)

# Pretty plot
mauc <- apply(array(unlist(auc), dim = c(nrow(auc[[1]]), ncol(auc[[1]]), length(auc))), c(1,2), mean)
pdf("auc.pdf", width=5, height=5)
matplot(log(lambda), mauc, type="l", lty=1, lwd=2, main="Breast cancer 10-year metastasis prognosis", ylab="AUC",xlab="log(lambda)")
legend("topright", legend=c("Gaussian", "SUQUAN"), col=seq(niter), lty=1, lwd=2)
dev.off()
