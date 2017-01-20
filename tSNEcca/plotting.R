N <- 1000 # Overal number of examples (train+test)
N_paired <- 500 # Number of training examples
MaxAngle <- 4*pi
MinRadius <-0.3
MaxRadius <- 8
NumNNs_X <- 20
NumNNs_Y <- 20
sx <- 0.5
sy <- 0.5
set.seed(8409)
## Generate data for views 1,2
t <- seq(0, MaxAngle, length.out = N)
r <- seq(MinRadius, MaxRadius, length.out = N) + 2*runif(N)
#### generate X, the noise can be added!
X <- cbind(r*cos(t+0*rnorm(N)*0.05), 
           r*sin(t+0*rnorm(N)*0.05)) 
X <- X + 0*matrix(rnorm(N*2), ncol = 2) 
#### generate Y, the noise can be added!
Y <- cbind(t+0*rnorm(N)*1, 
           2*rnorm(N))
Y <- Y + 0*cbind(rep(0, N), rnorm(N))
# ## Generate random data sets
# X <- apply(X, 2, sample)
# Y <- apply(Y, 2, sample)
## Training data
PairedIndices <- sample(1:N, N_paired)
# X <- read.table("X.dat", h=F, stringsAsFactors = F)
# Y <- read.table("Y.dat", h=F, stringsAsFactors = F)
# PairedIndices <- read.table("PairedIndices.dat", h=F, stringsAsFactors = F)
# 
# XX <- as.matrix(X); YY <- as.matrix(Y); PairedIndices <- as.numeric(PairedIndices)
## Test (or validation) data
UnpairedIndices <- setdiff(1:N,PairedIndices)
## plot data
par(mfrow = c(1, 2))
col_b2y <- colorRampPalette(c("blue", "cyan", "orange", "yellow"))
plot(X, pch = 19, col = col_b2y(N), 
     main = "View of X", xlab = "X.1", ylab = "X.2", xlim = range(X[,1]*1.1), ylim = range(X[,2]*1.1))
points(X[PairedIndices,], col = "navy")
legend("topright", c("Train", "Test"), pch = c(1, 19), col = "navy")
plot(Y, pch = 19, col = col_b2y(N), 
     main = "View of Y", xlab = "Y.1", ylab = "Y.2", xlim = range(Y[,1]*1.1), ylim = range(Y[,2]*1.3))
points(Y[PairedIndices,], col = "navy")
legend("topright", c("Train", "Test"), pch = c(1, 19), col = "navy")
par(mfrow = c(1,1))
