#####################
#
# Functions used to generate landscape characteristics
# Current version is minimal, with random location of 
# plots and random distribution of the environmental conditions
# Dominique Gravel
# June 2nd, 2016
# 
####################

#' Random XY coordinates
#' @param n number of patches
get_XY = function(N) cbind(runif(N),runif(N))

# Aggregation of XY coordinates
get_XY_agg = function(N, Nclusters, sd_xy) {

	Xclust = runif(Nclusters)
	Yclust = runif(Nclusters)

	X = rnorm(N, rep(Xclust,N/Nclusters), sd_xy)
	Y = rnorm(N, rep(Yclust,N/Nclusters), sd_xy)

	cbind(X,Y)
}

#' Random uniform environmental values
#' @param D number of environmental variables
#' @param N number of patches
#' @param SA type of autocorrelation. Default is none
get_E = function(D, N, SA = NULL, XY = NULL){
  if (is.null(SA)) out = matrix(runif(D*N), nr = N, nc = D)
  else{
    out = matrix(nr = N, nc = D)
    for (i in 1:D){
      Chol.matrix<-Chol.Spatial.Covariance(XY,range=SA[[1]],C0=2,C1=12)
      out[,i]<-as.vector(runif(N)%*%Chol.matrix)   
      out[,i] = out[,i] / max(out[,i])
    }
  }
  return(out) 
}
