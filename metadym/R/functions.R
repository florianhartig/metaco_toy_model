#####################
#
# Functions used to compute the colonization and extinction probabilities
# Names of functions correspond to the ones described in the main text
# These could be tweaked to change things such as the shape of the dispersal kernel
# or the response to the environment
# Dominique Gravel
# June 2nd, 2016
# 
####################

# Compute the sum of ecological interactions for every location and every species
sum_interactions = function (A, Y) t(A%*%t(Y))

# Compute the propagule pressure
I_f = function(Y, K, m) I = (1-m)*(K%*%Y)/(K%*%matrix(1,nr=N,nc=R)) + m

# Compute the local performance of propagules
S_f = function(E, u_c, s_c) {
	R = ncol(u_c)
	N = nrow(E)
	D = ncol(E)
	S = matrix(1, nr = N, nc = R)
	for(i in 1:D) S = S*exp(-(E[,i]-matrix(u_c[i,],nr=N,nc=R,byrow=TRUE))^2 / matrix(s_c[i,],nr=N,nc=R,byrow=TRUE)^2)
	return(S)
}

# Effect of ecological interactions on colonization probability
C_f = function(v, d_c, c_0, c_max) c_max*(1 +(1/c_0 - 1)*exp(-v*d_c))^-1

# Effect of the environment on the extinction
M_f = function(E, u_e, s_e) {
	R = ncol(u_e)
	N = nrow(E)
	D = ncol(E)
	M = matrix(1, nr = N, nc = R)
	for(i in 1:D) M = M*exp(-(E[,i]-matrix(u_e[i,],nr=N,nc=R,byrow=TRUE))^2 / matrix(s_e[i,],nr=N,nc=R,byrow=TRUE)^2)
	return(M)	
}

# Effect of ecological interactions on extinction
E_f = function(v, d_e, e_0, e_min) {

	e_min_mat = matrix(e_min, nr = N, nc = R, byrow=TRUE)

	e_min_mat+(1/(1-e_min_mat)+(1/(e_0-e_min_mat)-1/(1-e_min_mat))*exp(d_e*v))^-1

}
# Compute the connectivity matrix
get_K = function(XY, alpha) {
	N = nrow(XY)
	distMat = as.matrix(dist(XY, method = "euclidean", upper = T, diag = T))
	ConMat = exp(-1/alpha*distMat)
	diag(ConMat) = 0
	return(ConMat)
}


#' Compute the connectivity matrix
#' @param XY location of patches
#' @param alpha distance decay of dispersal
#' @param eps cutoff for sparse matrix
get_K_sparse = function(XY, alpha, eps = 0.001) {
  N = nrow(XY)
  distMat = as.matrix(dist(XY, method = "euclidean", upper = T, diag = T))
  ConMat = exp(-1/alpha*distMat)
  diag(ConMat) = 0
  
  nonZero = which(ConMat > eps, arr.ind = T)
  
  nonZeroValues = rep(NA, nrow(nonZero))
  
  for (i in 1:nrow(nonZero)){
    nonZeroValues[i] = ConMat[nonZero[i,1], nonZero[i,2]]
  }
  
  sparseConMat = Matrix::sparseMatrix(i = nonZero[,1], j = nonZero[,1], x = nonZeroValues)
  
  
  return(sparseConMat)
}

# XY = get_XY(100)
# alpha = 0.01
# eps = 0.001
# K = get_K(XY, 0.01)
# 
# K%*%XY


