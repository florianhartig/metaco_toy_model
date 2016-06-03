RDA <-function (y,X) 
{
  DistMatrix <- as.matrix(y)
  X <- as.matrix(X)
  n <- nrow(DistMatrix)
  p <- ncol(X)
  X<-X-mean(X)
  y<-y-mean(y)
  
  H<-X%*%solve(t(X)%*%X)%*%t(X)
  I<-diag(n)
  predicted <- H%*%y
  residuals <- y-predicted
  MS_regression<-sum(diag(t(predicted)%*%predicted))/p
  MS_residual<-sum(diag(t(residuals)%*%residuals))/(n-p)
  F<-MS_regression/MS_residual
  MS_Total=sum(diag(t(y)%*%y))/n;
  RsqAdj=1-MS_residual/MS_Total;
  
  result <- list(Rsq_Adj=RsqAdj,F_value=F,res_matrix=residuals,pred_matrix=predicted)
  return(result)
}

