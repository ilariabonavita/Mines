# tSNE-cca method
# try 
# description: Compute tSNE-correlation analysis between two data matrices
# 
wd.path <- '/Users/ilaria_bonavita/My_project/tSNEcca'
setwd(wd.path)

# dimension of the input features space
q <- 3
# dimension of the output transformed features
d <- 2
# Overal number of examples (train+test)
N <- 100
# Number of training examples
N_paired <- 100

training <- F

### Generate data 
MaxAngle <- 4*pi
MinRadius <-0.3
MaxRadius <- 8
sx <- 0.5
sy <- 0.5
set.seed(1244)
## Generate data for views 1,2
t <- seq(0, MaxAngle, length.out = N)
r <- seq(MinRadius, MaxRadius, length.out = N) + 2*runif(N)
#### generate X, the noise can be added!
set.seed(123)
#X <- cbind(rep(1:length(cos(t))),cos(t),rnorm(N))
#X <- cbind(r*sin(t+0*rnorm(N)*0.05),t,rnorm(N))
X <- cbind(r*cos(t+0*rnorm(N)*0.05),r*sin(t+0*rnorm(N)*0.05),rnorm(N))
#X <- X + 0*matrix(rnorm(N), ncol = q)
#### generate Y, the noise can be added!
#Y <- cbind(rep(1:length(sin(t))),sin(t),rnorm(N))
#Y <- cbind(r*cos(t+0*rnorm(N)*0.05),t,rnorm(N))
#Y <- cbind(r*sin(t+0*rnorm(N)*0.05),r*cos(t+0*rnorm(N)*0.05),rnorm(N))
Y <- cbind(t+0*rnorm(N)*1, 2*rnorm(N),rnorm(N))
#Y <- Y + 0*cbind(rep(0, N), rnorm(N))

if(anyNA(X)| anyNA(Y)){
  remove.na <- unlist(apply(cbind(X,Y), 1, anyNA))
  X <- X[!remove.na, , drop = F]
  Y <- Y[!remove.na, , drop = F]
  rm(remove.na) 
}

if (training){
## Training data
PairedIndices <- sample(1:N, N_paired)
## Test (or validation) data
UnpairedIndices <- setdiff(1:N,PairedIndices)
# Work with the train set
XTr <- X[PairedIndices,]
YTr <- Y[PairedIndices,]
}

XTr <- X
YTr <- Y

#### Define a matrix of the distance between Xi and Xj (Yi and Yj) for all the possible combination of elements in X (Y)
XTr.pairs <- combn(seq(1:nrow(XTr)),2)
YTr.pairs <- combn(seq(1:nrow(YTr)),2)

# For simplicity, we use only the distances from now on
## Each row correspond to a pair (i,j), columns are the coordinates of the distance xi-xj
XTr.dist <- XTr[XTr.pairs[1,],] - XTr[XTr.pairs[2,],]
YTr.dist <- YTr[YTr.pairs[1,],] - YTr[YTr.pairs[2,],]

#=================================================================
#############  Batch gradient descent
source('./batch.grad.desc.fun.R')
batch.grad.desc.fun(XTr.dist, YTr.dist, N, dim.out=2,
                                gamma.list= c(0.06, 0.05, 0.02, 0.005),
                                nexper = 1,
                                maxiter= 500,
                                perf.cca = T,
                                wd.path= '/Users/ilaria_bonavita/My_project/tSNEcca'
                                )
batch.grad.desc.fun(XTr.dist, YTr.dist, N, dim.out=2,
                    gamma.list= c(0.06,0.05),
                    nexper = 1,
                    maxiter= 200,
                    perf.cca = F,
                    wd.path= '/Users/ilaria_bonavita/My_project/tSNEcca',
                    exp.name='cos_sin_const'
                    )

###==============================================================================
### EXAMPLE with plots

# Retrieve the projection matrices
mat.path <- '/Users/ilaria_bonavita/My_project/tSNEcca/N100_in3_out2_rnd/exp_1gamma0.05/'
Wout <- readRDS(paste0(mat.path,'Wout.RDS'))
Zout <- readRDS(paste0(mat.path,'Zout.RDS'))

source('./tsne2Dplot.fun.R')
tsne2Dplot.fun(XTr,YTr,Wout,Zout,'/Users/ilaria_bonavita/My_project/tSNEcca/N100_in3_out2_cos_sin_const/exp_1gamma0.05')
tsne2Dplot.fun(XTr,YTr,Wout.nocca,Zout.nocca,'/Users/ilaria_bonavita/My_project/tSNEcca/N100_in3_out2_rnd/exp_1gamma0.05')


# Random W and Z
#set.seed(151)

#Wrand <- matrix(rnorm(N*ncol(XTr.dist)),nrow = d, ncol = ncol(XTr.dist))
#Wrand <- t(apply(Wrand,1, function(x){ x/sum(x)}))
#Zrand <- matrix(rnorm(N*ncol(YTr.dist)),nrow = d, ncol = ncol(YTr.dist))
#Zrand <- t(apply(Zrand,1, function(x){ x/sum(x)}))
 
#XWrand <- Wrand%*%t(X[PairedIndices,])
#YZrand <- Zrand%*%t(Y[PairedIndices,])          

#pdf(paste0(mat.path,'random_projections.pdf'))
#par(mfrow = c(1, 2))
#col_b2y <- colorRampPalette(c("blue", "cyan", "orange", "yellow"))
#plot(t(XWrand), pch = 19, col = col_b2y(N), 
#     main = "View of X", xlab = "X.1", ylab = "X.2", xlim = range(XWrand[1,]*1.1), ylim = range(XWrand[2,]*1.1))
#plot(t(YZrand), pch = 19, col = col_b2y(N), 
#     main = "View of Y", xlab = "Y.1", ylab = "Y.2", xlim = range(YZrand[1,]*1.1), ylim = range(YZrand[2,]*1.3))
# dev.off()
 
 
# Perfomr NCCA
require(FNN)
require(Matrix)
require(irlba)
source('./ncca.R')
ncca_res <- ncca(X,Y, d = 2, hx = 0.75, hy = 0.75, nx=5, ny=5)

X_proj_paired <- ncca_res$X_new
Y_proj_paired <- ncca_res$Y_new

x11()
par(mfrow = c(1, 2))
col_b2y <- colorRampPalette(c("blue", "cyan", "orange", "yellow"))
plot(X_proj_paired, pch = 19, col = col_b2y(N), 
     main = "View of X", xlab = "X.1", ylab = "X.2", xlim = range(X_proj_paired[,1]*1.1), ylim = range(X_proj_paired[,2]*1.1))
plot(Y_proj_paired, pch = 19, col = col_b2y(N), 
     main = "View of Y", xlab = "Y.1", ylab = "Y.2", xlim = range(Y_proj_paired[,1]*1.1), ylim = range(Y_proj_paired[,2]*1.3))

x11()
par(mfrow = c(1, 2))
col_b2y <- colorRampPalette(c("blue", "cyan", "orange", "yellow"))
plot(X_proj_paired[,1], Y_proj_paired[,1] ,pch = 19, col = col_b2y(N), 
     main = "View of X1 vs Y1", xlab = "X.1", ylab = "Y.1", xlim = range(X_proj_paired[,1]*1.1), ylim = range(Y_proj_paired[,1]*1.1))
plot(X_proj_paired[,2],Y_proj_paired[,2], pch = 19, col = col_b2y(N), 
     main = "View of X2 vs Y2", xlab = "X.2", ylab = "Y.2", xlim = range(X_proj_paired[,2]*1.1), ylim = range(Y_proj_paired[,2]*1.3))


# Perform KCCA
library(kernlab)
library(geigen)
source('./kerncca.R')

kerncca_res <- kerncca(X,Y)

x11()
par(mfrow = c(1, 2))
col_b2y <- colorRampPalette(c("blue", "cyan", "orange", "yellow"))
plot(kerncca_res$latent1, pch = 19, col = col_b2y(N), 
     main = "View of X", xlab = "X.1", ylab = "X.2", xlim = range(kerncca_res$latent1[,1]*1.1), ylim = range(kerncca_res$latent1[,2]*1.1))
plot(kerncca_res$latent2, pch = 19, col = col_b2y(N), 
     main = "View of Y", xlab = "Y.1", ylab = "Y.2", xlim = range(kerncca_res$latent2[,1]*1.1), ylim = range(kerncca_res$latent2[,2]*1.3))

x11()
par(mfrow = c(1, 2))
col_b2y <- colorRampPalette(c("blue", "cyan", "orange", "yellow"))
plot(kerncca_res$latent1[,1], kerncca_res$latent2[,1] ,pch = 19, col = col_b2y(N), 
     main = "View of X1 vs Y1", xlab = "X.1", ylab = "Y.1", xlim = range(kerncca_res$latent1[,1]*1.1), ylim = range(kerncca_res$latent2[,1]*1.1))
plot(kerncca_res$latent1[,2],kerncca_res$latent2[,2], pch = 19, col = col_b2y(N), 
     main = "View of X2 vs Y2", xlab = "X.2", ylab = "Y.2", xlim = range(kerncca_res$latent1[,2]*1.1), ylim = range(kerncca_res$latent2[,2]*1.3))


# Perform CCA
cca_res  <- cancor(X,Y)

XWcca <- X %*% cca_res$xcoef
YZcca <- Y %*% cca_res$ycoef

x11()
par(mfrow = c(1, 2))
col_b2y <- colorRampPalette(c("blue", "cyan", "orange", "yellow"))
plot(XWcca, pch = 19, col = col_b2y(N), 
     main = "View of X", xlab = "X.1", ylab = "X.2", xlim = range(XWcca[,1]*1.1), ylim = range(XWcca[,2]*1.1))
plot(YZcca, pch = 19, col = col_b2y(N), 
     main = "View of Y", xlab = "Y.1", ylab = "Y.2", xlim = range(YZcca[,1]*1.1), ylim = range(YZcca[,2]*1.3))

x11()
par(mfrow = c(1, 2))
col_b2y <- colorRampPalette(c("blue", "cyan", "orange", "yellow"))
plot(XWcca[,1], YZcca[,1] ,pch = 19, col = col_b2y(N), 
     main = "View of X1 vs Y1", xlab = "X.1", ylab = "Y.1", xlim = range(XWcca[,1]*1.1), ylim = range(YZcca[,1]*1.1))
plot(XWcca[,2],YZcca[,2], pch = 19, col = col_b2y(N), 
     main = "View of X2 vs Y2", xlab = "X.2", ylab = "Y.2", xlim = range(XWcca[,2]*1.1), ylim = range(YZcca[,2]*1.3))

