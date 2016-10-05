#--- library loading
library(methods)
library(kernlab) # for kernel svm
library(pcaPP) # for fast kendall tau
library(caret) # for data split
library(parallel) # for parallel cross-validation
library(kernrank)

setwd('/Users/ilaria_bonavita/Mines/mines_yunlong')
source('/Users/ilaria_bonavita/Mines/mines_yunlong/func.R')

sink('log.txt',append = T)

#=============================== Dataset
#--- load data
data.path <- './data/arloth.txt'
dataset <- read.table(data.path,h=T)

#--- give the right format  (NOTE no NA, NaN are present)
mdd <- list(xtrain = as.matrix(dataset[,-1]), ytrain= dataset[,1], xtest=NULL, ytest=NULL)
mdd$ytrain <- as.factor(ifelse(mdd$ytrain==1,'case','control'))

#--- balance the dataset 
ind1 <- which(mdd$ytrain %in% 'case')
ind0 <- which(mdd$ytrain %in% 'control')
sampsize <- min(length(ind1),length(ind0))
sampind1 <- sample(ind1, sampsize)
sampind0 <- sample(ind0, sampsize)
sampind <- c(sampind0,sampind1)
mdd.bal <- list(xtrain = mdd$xtrain[sampind,], ytrain= mdd$ytrain[sampind],xtest=NULL, ytest=NULL)

datasetlist <- c('mdd','mdd.bal')

#=============================== Model performance comparison
#--- set list of C parameters 
Cpara_list <- 10^(-2:3)
names(Cpara_list) <- paste('C',1:length(Cpara_list),sep='')

#--- categorize models for ease of training
modelsConly <- c("SVMlinearALL", "SVMkdtALL", "SVMpolynomialALL", "SVMrbf") # plus KFD coded within each


#--- training loop
modelsKFD <- sub("SVM", "KFD", modelsConly)

  ptm <- proc.time()
  
  #--- crossval datasets
  res_crossval <- mclapply(datasetlist, function(datasetlist){
    xtr <- get(datasetlist)$xtrain; ytr <- get(datasetlist)$ytrain
    xtst <- get(datasetlist)$xtest; ytst <- get(datasetlist)$ytest
    if(!is.null(xtst) || !is.null(ytst)) stop(paste('dataset error',datasetlist,sep=':'))
    
    set.seed(1226)
    outterFoldIndices <- createMultiFolds(1:nrow(xtr), k=5, times=10)
    sig <- sigest(xtr,scaled=F)['50%']
    
    res <- lapply(outterFoldIndices, function(outterFold){
      return(perfClassification(NULL, datasetlist, xtr[outterFold,,drop=F], ytr[outterFold], xtr[-outterFold,,drop=F], ytr[-outterFold], 
                                Cpara_list, NULL, modelsConly, NULL, NULL, 
                                nfolds = 5, nrepeats = 1, seed =(206), sigma=sig))
    })
    return(res)
  }, mc.cores = 3,mc.preschedule = FALSE)
  names(res_crossval) <- datasetlist
  
  #==============================
  table_acc <- matrix(-100, 
                      nrow = length(datasetlist), ncol = length(c(modelsConly,modelsKFD)),
                      dimnames = list(datasetlist, c(modelsConly,modelsKFD)))
  
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
write.table(table_acc,'table_acc.txt',  qu=F, col.names = T)

sink()


