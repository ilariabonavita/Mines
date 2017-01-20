
batch.grad.desc.fun <- function(X.dist,
                                Y.dist,
                                Nsample,
                                dim.out=2,
                                gamma.list= c(0.05, 0.01, 0.005, 0.001),
                                nexper = 1,
                                maxiter= 200,
                                perf.cca=F,
                                wd.path= '/Users/ilaria_bonavita/My_project/tSNEcca',
                                exp.name=NULL
                                )
{ 
  exp.path <- paste0(wd.path,'/N',Nsample,'_in',ncol(X.dist),'_out',dim.out)
  if(!is.null(exp.name))
    exp.path <- paste0(exp.path,'_',exp.name)
  
  dir.create(exp.path)
  
  for (k in 1:nexper){
    
    # Inizialize W,Z for each experiment
    if (perf.cca){
    # iniziatlize from cca  
      library('CCA')
      print('CCA performed on the distances matrices to initialize W and Z')
      cca.out <- cc(XTr.dist,YTr.dist)
      W <- cca.out$xcoef
      Z <- cca.out$ycoef
    }  else { 
    # random initialization
    set.seed(142+k)
    W <- matrix(rnorm(Nsample*ncol(X.dist)),nrow = dim.out, ncol = ncol(X.dist))
    Z <- matrix(rnorm(Nsample*ncol(Y.dist)),nrow = dim.out, ncol = ncol(Y.dist))
    }
    # normalize projection matrices  
    W <- t(apply(W,1, function(x){ x/sum(x)}))
    Z <- t(apply(Z,1, function(x){ x/sum(x)}))
    
    W.base <- W
    Z.base <- Z
    Wold <- W.base
    Zold <- Z.base
    
    for (gamma in gamma.list){
      fold.path <- paste0(exp.path,'/exp_',k,'gamma',gamma)
      dir.create(fold.path)
      # Create a vector where to save the difference C.new - C.old at each iteration 
      diff <- rep(0, maxiter)
      
      for(i in c(1: maxiter)){
        source(paste0(wd.path,'/distr_cost_fun.R'))
        source(paste0(wd.path,'/gradient_functions.R'))
        
        allnum.P <- allnumerator.vec(X.dist,Wold)
        allnum.Q <- allnumerator.vec(Y.dist,Zold)
        
        distr.P <- distr.vec(X.dist,Wold,allnum.P)
        distr.Q <- distr.vec(Y.dist,Zold,allnum.Q)
        
        C.dW <- C.deriv(X.dist,Wold,allnum.P,distr.P,distr.Q)
        C.dZ <- C.deriv(Y.dist,Zold,allnum.Q,distr.Q,distr.P)
        
        W <- Wold - gamma*C.dW
        Z <- Zold - gamma*C.dZ
        
        #print(W)
        #print(Z)
        
        #--- compute difference Cold Cnew
        allnum.P.new <- allnumerator.vec(X.dist,W)
        allnum.Q.new <- allnumerator.vec(Y.dist,Z)
        distr.P.new <- distr.vec(X.dist,W,allnum.P.new)
        distr.Q.new <- distr.vec(Y.dist,Z,allnum.Q.new)
        
        diff[i] <- C.fun(distr.P,distr.Q) - C.fun(distr.P.new,distr.Q.new)
       # diff[i] <-  C.fun(distr.P.new,distr.Q.new)
        #print(diff)
        
        #-- update Wold
        Wold <- W
        Zold <- Z
      }
      pdf(paste0(fold.path,'/convergence.pdf'))
      plot(1:maxiter,diff, xlab = 'iteration', ylab = 'delta cost function',main = paste0("Batch Gradient Descendent-gamma=",gamma))
      dev.off()
      C.min <-  C.fun(distr.P.new,distr.Q.new)
      saveRDS(C.min,paste0(fold.path,'/Cmin.RDS'))
      saveRDS(Wold,paste0(fold.path,'/Wout.RDS'))
      saveRDS(Zold,paste0(fold.path,'/Zout.RDS'))
      saveRDS(diff,paste0(fold.path,'/diff.RDS'))
    }
  }
 
}