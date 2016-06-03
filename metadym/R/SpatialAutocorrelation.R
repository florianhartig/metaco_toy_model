Chol.Spatial.Covariance<-function(Coordinates,range=15,C0=2,C1=12){
  u<-as.matrix(dist(Coordinates))
  n<-dim(Coordinates)[1]
  R<-matrix(0,n,n)  
  for (i in 1:n){
    for (j in 1:n){
      if ((u[i,j]<=range)&(i<=j)){
        R[i,j]<-1-1.5*(u[i,j]/range)+0.5*((u[i,j]/range)^3)
        R[j,i]<-R[i,j]
      }
    }
  }
  Sigma<-C0*diag(dim(Coordinates)[1])+C1*R
  Chol<-chol(Sigma)
  return(Chol)
}

# example
#rm(list = ls())
