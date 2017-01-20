# Script for checking the gradients 

gradcheck <- function(W,Z,X.dist,Y.dist,eps,i,j){
  source("/Users/ilaria_bonavita/My_project/tSNEcca/distr_cost_fun.R")
  source("/Users/ilaria_bonavita/My_project/tSNEcca/gradient_functions.R")
  
  W.test <- W
  W.plus <- W.test
  W.minus <- W.test
  W.plus[i,j] <- W.test[i,j] + eps
  W.minus[i,j] <- W.test[i,j] - eps
  
  allnum.P<- allnumerator.vec(X.dist,W.test)
  allnum.P.plus <- allnumerator.vec(X.dist,W.plus)
  allnum.P.minus <- allnumerator.vec(X.dist,W.minus)
  distr.P<- distr.vec(X.dist,W.plus,allnum.P)
  distr.P.plus <- distr.vec(X.dist,W.plus,allnum.P.plus)
  distr.P.minus <- distr.vec(X.dist,W.minus,allnum.P.minus)
  allnum.Q <- allnumerator.vec(Y.dist,Z)
  distr.Q <- distr.vec(Y.dist,Z,allnum.Q)
  
  dP.approx <- (distr.P.plus - distr.P.minus)/(2*eps)
  dP.exact <- deriv.3dmat(X.dist,W.test,allnum.P,distr.P)
  
  C.deriv.approx <- (C.fun(distr.P.plus,distr.Q) - C.fun(distr.P.minus,distr.Q))/(2*eps)
  print(C.deriv.approx)
  C.deriv <- C.deriv(X.dist,W.test,allnum.P,distr.P,distr.Q)
  print(C.deriv)
  
  dCout <- ((C.deriv.approx - C.deriv[i,j])/max(C.deriv.approx,C.deriv[i,j]))<1e-4
  
  return(dCout)
}



