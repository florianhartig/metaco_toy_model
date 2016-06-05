#' Main function performing the spatially explicit model
#' @note  All of the different steps are functions, provided in other scripts
#' @author  Dominique Gravel
runModel = function(pars, nsteps) {
  
  
  # Landscape characteristics
  XY = metadym::get_XY(N)
  E = metadym::get_E(D,N, SA = pars$SA, XY = XY)
  
  # Initial conditions
  Y0 = matrix(0, nr = N, nc = R)
  rand = matrix(runif(N*R), nr = N, nc = R)
  Y0[rand < 0.5] = 1
  

  with(pars, {

	# Simulation parameters
	N = nrow(XY)
	R = ncol(u_c)

	# Initial conditions 
	Y = Y0 

	# Compute the connectivity matrix
	K = metadym::get_K(XY, alpha)

	# Compute the local performance
	S = S_f(E, u_c, s_c)
	M = M_f(E, u_e, s_e) 
  e0_E = (1 - matrix(e_0,nr=N,nc=R,byrow=TRUE))*(1-M) + matrix(e_0,nr=N,nc=R,byrow=TRUE)

  # Store the results
  
	RES = array(dim=c(nsteps,N,R))

  	# Main loop
  	for(time in 1:nsteps) {

  		# Matrix to store changes
  		delta = matrix(0, nr = N, nc = R)

  		### Colonization ###
  		# Compute elements of the colonization probability
  		v = sum_interactions(A, Y)
  		I = I_f(Y, K, m)
  		C = C_f(v, d_c, c_0, c_max)

  		# Colonization prob
  		P_col = I*S*C

  		# Perform the test
  		rand = matrix(runif(N*R), nr = N, nc = R)
  		delta[Y == 0 & rand < P_col] = 1

  		### Extinction ###
  		# Compute the extinction probability
  		P_ext = E_f(v, d_e, e_0 = e0_E, e_min)

  		# Perform the test
  		rand = matrix(runif(N*R), nr = N, nc = R)
  		delta[Y == 1 & rand < P_ext] = - 1	

  		### Apply changes ###
  		Y = Y + delta

  		### Record results ###
  		RES[time,,] = Y
  	 
 #    cat("Step = ", time, '\n')
    } # End of loop 
	
	
	
  return(list(occupancy = RES, XY = XY, E = E, pars = pars))
  })
}






