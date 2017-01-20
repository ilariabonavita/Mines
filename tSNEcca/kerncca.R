#' Kernel Canonical Correlation Analysis
#'
#' @description Compute the kernel canonical correlations between two data matrices
#'
#' @importFrom kernlab inchol()
#' @importFrom geigen geigen()
#'
#' @param x1 The data set containing n samples and p features
#' @param x2 The data set containing n samples and q features
#' @param ktyp1 The kernel for the data set x1 (The default value is "rbfdot".)
#' @param ktyp2 The kernel for the data set x2 (The default value is "rbfdot".)
#' @param kpar1 The list of kernel parameters for kernel 1
#' @param kpar2 The list of kernel parameters for kernel 1
#' @param reg The regularization parameter (The default is 1e-5.)
#' @param ...
#' @details The incomplete Cholesky decomposition is calculated by inchol() from kernlab. The details of choosing the kernel function and setting the kernel parameters is seen in \url{http://finzi.psych.upenn.edu/library/kernlab/html/inchol.html}
#'
#' @return cor The kernel canconical correlation between x1 and x2, which is the first canonical correlation between projections of x1 and x2
#' @return latent1 The latent variable, which is the nonlinear projection of x1 in the kernel feature space
#' @return latent2 The latent variable, which is the nonlinear projection of x2 in the kernel feature space
#'
#' @export
#' @author Meiwen Jia, \email{meiwen_jia@psych.mpg.de}
#' @references Vaerenbergh, S. V. (2010). Kernel methods for nonlinear identification, equalization and separation of signals. Universidad de Cantabria.
#'
#' @examples
#' set.seed(111)
#' N <- 1000
#' x <- runif(N) # latent signal
#' r1 <- runif(N);  r2 <- runif(N)# random (helper) variables
#' x1 <- cbind(tan(r1-x)+0.1*r1, r1+3*x-1/10*(sin(3*x)))
#' x2 <- cbind(x - 2*(1-exp(-r2))/(1+exp(-r2)), r2*x, tan(r2+x))
#' cat("The correlation is ", kerncca(x1, x2)$cor, "\n")
#'
kerncca <- function(x1,
                    x2,
                    ktyp1 = "rbfdot",
                    ktyp2 = "rbfdot",
                    kpar1 = list(sigma=0.01),
                    kpar2 = list(sigma=0.01),
                    reg = 1e-5, # regularization
                    nlatent =2 ,
                    ...){
  # An incomplete cholesky decomposition calculates Z where K= ZZ' K being the kernel matrix.
  G1 <- kernlab::inchol(x1, kernel = ktyp1, kpar = kpar1)
                        #, ...)
  G2 <- kernlab::inchol(x2, kernel = ktyp2, kpar = kpar2)
  #, ...)
  # do column centering
  # important! ignoring will cause problem
  G1 <- scale(G1, center = T, scale = F)
  G2 <- scale(G2, center = T, scale = F)
  #   # option1 faster!
  #   # ones and zeros
  #   N1 <- ncol(G1)
  #   N2 <- ncol(G2)
  #   Z11 <- matrix(0, N1, N1)
  #   Z22 <- matrix(0, N2, N2)
  #   Z12 <- matrix(0, N1, N2)
  #   I11 <- diag(N1)
  #   I22 <- diag(N2)
  #
  #   # simplified Hardoon
  #   R <- rbind(cbind(Z11, t(G1)%*%G2),
  #              cbind(t(G2)%*%G1, Z22))
  #   D <-rbind(cbind(t(G1)%*%G1+reg*I11,Z12),
  #             cbind(t(Z12), t(G2)%*%G2+reg*I22));
  #
  #   # solve generalized eigenvalue problem
  #   eig.res <- geigen::geigen(R, D)
  #   vecs <- eig.res$vectors
  #   vals <- eig.res$values
  #   alpha <- vecs[, which.max(vals), drop = F]
  #   alpha <- alpha/norm(alpha, "1")
  #   # the coefficient is the maximum eigenvalues
  #   beta <- max(vals)
  #
  #   # expansion coefficients
  #   alpha1 <- matrix(alpha[1:N1], ncol=1, byrow = F);
  #   alpha2 <- matrix(alpha[N1+(1:N2)], ncol=1, byrow = F);
  #
  #   # estimates of latent variable
  #   y1 <- G1%*%alpha1
  #   y2 <- G2%*%alpha2
  #
  #
  #   return(list(cor = sort(vals, decreasing = T), latent1 = y1, latent2 = y2, G1 = G1, G2 = G2))
  # option2
  G1G2_cca <- cancor(G1, G2)
  alpha <- matrix(c(G1G2_cca$xcoef[,1:nlatent], G1G2_cca$ycoef[,1:nlatent]),ncol = nlatent)
                  #, ncol = 1)
  alpha <- alpha/norm(alpha, "1")
  y1 <- G1%*%alpha[1:ncol(G1),1:nlatent]
  y2 <- G2%*%alpha[ncol(G1)+(1:ncol(G2)),1:nlatent]
  return(list(cor = G1G2_cca$cor, latent1 = y1, latent2 = y2, G1 = G1, G2 = G2))
}
