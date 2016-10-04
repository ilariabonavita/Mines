# utiles
library(kernlab) # for kernel svm
library(pcaPP) # for fast kendall tau
library(caret) # for data split
library(parallel) # for parallel cross-validation
library(kernrank)

setwd('/Users/ilaria_bonavita/Documents/Paris/R/kendallkernel_demo/geneexpr')
source("func.R")

#=============================== Dataset
#--- load data
data.path <- './data/arloth.txt'
dataset <- read.table(data.path,h=T)

#--- give the right format  (NOTE no NA, NaN are present)
mdd <- list(xtrain = as.matrix(dataset[,-1]), ytrain= dataset[,1], xtest=NULL, ytest=NULL)
mdd$ytrain <- as.factor(ifelse(mdd$ytrain==1,'case','control'))

#--- balance the dataset 
if (balance==T) {
  ind1 <- which(mdd$ytrain %in% 'case')
  ind0 <- which(mdd$ytrain %in% 'control')
  sampsize <- min(length(ind1),length(ind0))
  sampind1 <- sample(ind1, sampsize)
  sampind0 <- sample(ind0, sampsize)
  sampind <- c(sampind0,sampind1)
  mdd.bal <- list(xtrain = mdd$xtrain[sampind,], ytrain= mdd$ytrain[sampind],xtest=NULL, ytest=NULL)
}

datasetlist <- c('mdd','mdd.bal')

#=============================== Model performance comparison
#--- set list of C parameters for SVM-based models
Cpara_list <- 10^(-2:3)
names(Cpara_list) <- paste('C',1:length(Cpara_list),sep='')

#--- categorize models for ease of training
#modelsNOpara <- c("APMV")
modelsConly <- c("SVMlinearALL", "SVMkdtALL", "SVMpolynomialALL", "SVMrbf") # plus KFD coded within each
#modelsTSPrelated <- c("TSP", "kTSP", "SVMlinearTOP", "SVMkdtTOP", "SVMpolynomialTOP")
modelslist <- 'modelsConly'
# set list of #genepairs for corresponding models
max_nodes = 5000; npairs_out = 30;
npairs_list <- floor(exp(seq(0,1,length.out=npairs_out)*log(max_nodes)))
evenidx <- npairs_list %% 2 == 0
npairs_list[evenidx] <- npairs_list[evenidx] - 1 # keep odd numbers only
npairs_list <- unique(npairs_list)
names(npairs_list) <- paste('k',1:length(npairs_list),sep='')

#--- crossval datasets
res_crossval <- mclapply(datasetlist, function(datasetlist){
  xtr <- get(datasetlist)$xtrain; ytr <- get(datasetlist)$ytrain
  xtst <- get(datasetlist)$xtest; ytst <- get(datasetlist)$ytest
  if(!is.null(xtst) || !is.null(ytst)) stop(paste('dataset error',datasetlist,sep=':'))
  
  set.seed(1226)
  outterFoldIndices <- createMultiFolds(1:nrow(xtr), k=2, times=1)
  sig <- sigest(xtr,scaled=F)['50%']
  
  res <- lapply(outterFoldIndices, function(outterFold){
    return(perfClassification(NULL, datasetlist, xtr[outterFold,,drop=F], ytr[outterFold], xtr[-outterFold,,drop=F], ytr[-outterFold], 
                              Cpara_list, npairs_list, modelsConly, NULL, NULL, 
                              nfolds = 3, nrepeats = 1, seed = 206, sigma=sig))
  })
  return(res)
}, mc.cores = 3,mc.preschedule = FALSE)
names(res_crossval) <- datasetlist

#==============================
modelsKFD <- sub("SVM", "KFD", modelsConly)
table_acc <- matrix(-100, 
                    nrow = length(datasetlist), ncol = length(c(modelslist,modelsKFD)),
                    dimnames = list(datasetlist, c(modelslist,modelsKFD)))

for (prefixname in datasetlist) {
  for (modelname in modelsConly) {
        s <- mean(sapply(res_crossval[[prefixname]], function(res){
        idx <- which.max(res[[modelname]]$cvacc)
        return(res[[modelname]]$acc[idx])
      }))
        
      table_acc[prefixname,modelname] <- round(100*s,2)
      s_kfd <- mean(sapply(res_crossval[[prefixname]], function(res){
          return(res[[modelname]]$acc_kfd)
        }))
      
        table_acc[prefixname,sub("SVM", "KFD", modelname)] <- round(100*s_kfd,2)

    }
  }

rownames(table_acc) <- nameafter[match(rownames(table_acc), namebefore)] # re-name
table_acc <- table_acc[order(rownames(table_acc)), ] # re-order
table_acc <- rbind(AVERAGE = round(colMeans(table_acc), 2), table_acc) # add AVERAGE scores over all datasets
table_acc <- table_acc[ ,order(table_acc["AVERAGE",],decreasing = TRUE)] # re-order
# show score table
t(table_acc)
# show boxplot
par(mar = c(10, 5, 1, 1) + 0.1, font.lab = 2, font.axis = 2, font.main = 2, cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5)
boxplot(table_acc[-1, ]/100, las = 2, ylab = 'cacc', col='royalblue2')

nConly <- length(modelsConly) # number of SVM models (KFD implemented within each)
for (prefixname in datasetlist) {
    mainname <- prefixname
    res <- res_crossval[[prefixname]]
    s <- lapply(res, function(resfold){
      sapply(resfold[modelsConly],function(u){c(u$acc_kfd,u$acc)})
    })
    dm <- c(dim(s[[1]]),length(s)); dn <- dimnames(s[[1]]) # save info of dim and dimnames
    s <- unlist(s); dim(s) <- dm # reform to an 3d array
    s <- apply(s, c(1,2), mean); dimnames(s) <- dn # average over cv-folds
  
  # set y range for plot
  plotrange <- c(floor(10*min(s,na.rm=T))/10,ceiling(10*max(s,na.rm=T))/10)
  
  # plotting
  par(font.lab = 2, font.axis = 2, font.main = 2, font = 2, cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5)
  plot(Cpara_list, rep(-100,length(Cpara_list)), main=mainname, 
       xlab="C parameter", ylab="acc", 
       ylim=plotrange,type='l',lwd=2,lty=1,col=1,log='x')
  for (col in seq(nConly)) {
    modelname <- modelsConly[col]
    lines(Cpara_list,s[-1,modelname],type='l',lwd=2,lty=1,col=col) # *** see notes above
    points(Cpara_list,rep(s[1,modelname],length(Cpara_list)),type='b',lty=5,lwd=1,pch=col,col=col,cex=1)
   
  }
  molist <- c(modelsConly,sub("SVM","KFD",modelsConly))
  molist[grep('rbf',molist)] <- paste(molist[grep('rbf',molist)],'ALL',sep='')
  legend("bottomright",legend=molist,pch=c(rep(NA,nConly),seq(nConly)),
         col=c(seq(nConly),seq(nConly)),lty=c(rep(1,nConly),rep(5,nConly)),
         lwd=c(rep(2,nConly),rep(1,nConly)),cex=1.25)
  grid(ny=16)
}
