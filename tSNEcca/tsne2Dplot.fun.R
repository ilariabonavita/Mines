tsne2Dplot.fun <- function(X,Y,W,Z,plot.path){
  
  N <- nrow(X)
  # Plot original data
  plot.basename <- paste0(plot.path,'/N',nrow(X),'_',nrow(W),'in_',ncol(W),'out')
  pdf(paste0(plot.basename,'_original.pdf'))
  par(mfrow = c(1, 2))
  col_b2y <- colorRampPalette(c("blue", "cyan", "orange", "yellow"))
  plot(X, pch = 19, col = col_b2y(N), 
       main = "View of X", xlab = "X.1", ylab = "X.2", xlim = range(X[,1]*1.1), ylim = range(X[,2]*1.1))
 # points(X, col = "navy")
  #legend("topright", c("Train", "Test"), pch = c(1, 19), col = "navy")
  plot(Y, pch = 19, col = col_b2y(N), 
       main = "View of Y", xlab = "Y.1", ylab = "Y.2", xlim = range(Y[,1]*1.1), ylim = range(Y[,2]*1.3))
  #points(Y, col = "navy")
  #legend("topright", c("Train", "Test"), pch = c(1, 19), col = "navy")
  dev.off()
  
  
  pdf(paste0(plot.basename,'_coordinates.pdf'))
  par(mfrow = c(1, 2))
  col_b2y <- colorRampPalette(c("blue", "cyan", "orange", "yellow"))
  plot(X[,1], Y[,1] ,pch = 19, col = col_b2y(N), 
       main = "X1 vs Y1", xlab = "X.1", ylab = "Y.1", xlim = range(X[,1]*1.1), ylim = range(Y[,1]*1.1))
  plot(X[,2], Y[,2] ,pch = 19, col = col_b2y(N), 
       main = "X2 vs Y2", xlab = "X.2", ylab = "Y.2", xlim = range(X[,2]*1.1), ylim = range(Y[,2]*1.1))
  dev.off()
  
  ### Plot transformed data
  XW <- W%*%t(X)
  YZ <- Z%*%t(Y) 
  
  # This plots the two views in the projected learned space - work only if outer space >1 
  x11()
  pdf(paste0(plot.basename,'_transformed.pdf'))
  par(mfrow = c(1, 2))
  col_b2y <- colorRampPalette(c("blue", "cyan", "orange", "yellow"))
  plot(t(XW), pch = 19, col = col_b2y(N), 
       main = "View of X", xlab = "X.1", ylab = "X.2", xlim = range(XW[1,]*1.1), ylim = range(XW[2,]*1.1))
  plot(t(YZ), pch = 19, col = col_b2y(N), 
       main = "View of Y", xlab = "Y.1", ylab = "Y.2", xlim = range(YZ[1,]*1.1), ylim = range(YZ[2,]*1.3))
  dev.off()
  
  pdf(paste0(plot.basename,'_coordinates_transformed.pdf'))
  par(mfrow = c(1, 2))
  col_b2y <- colorRampPalette(c("blue", "cyan", "orange", "yellow"))
  plot(XW[1,], YZ[1,] ,pch = 19, col = col_b2y(N), 
       main = "View of X1 vs Y1", xlab = "X.1", ylab = "Y.1", xlim = range(XW[1,]*1.1), ylim = range(YZ[1,]*1.1))
  plot(XW[2,],YZ[2,], pch = 19, col = col_b2y(N), 
       main = "View of X2 vs Y2", xlab = "X.2", ylab = "Y.2", xlim = range(XW[2,]*1.1), ylim = range(YZ[2,]*1.3))
  
  dev.off()
}