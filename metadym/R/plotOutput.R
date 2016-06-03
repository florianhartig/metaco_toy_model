#' Plots an overview of the model output
#' @param out model output
plotOutput <- function(out){
  
  par(mfrow = c(2,2))
  
  
  plot(out$XY[,1], out$XY[,2], cex = 2*out$E[,1], main = "Environment")
  
  div = computeAlphaGamma(out)
  alpha_div = div[,1]
  gamma_div = div[,2]
  
  #dev.new(width = 8, height = 6)
  par(mar = c(5,6,2,1))
  plot(c(1:nsteps), alpha_div, type = "l", xlab = "Time", ylab = "Diversity", 
       cex.lab = 1.25, cex.axis = 1.25)
  title(main = "Local diversity")
  
  #dev.new(width = 8, height = 6)
  par(mar = c(5,6,2,1))
  plot(c(1:nsteps), gamma_div, type = "l", xlab = "Time", ylab = "Diversity", 
       cex.lab = 1.25, cex.axis = 1.25)
  title(main = "Regional diversity")
  
  
  plot(out$colonization, type = "l", main = "Colonization")
  
}

computeAlphaGamma <- function(out){
  # Compute the number of species per time step
  nsteps = length(out$occupancy)
  alpha_div = numeric(nsteps)
  gamma_div = numeric(nsteps)
  
  for(i in 1:nsteps) {
    alpha_div[i] = mean(apply(out$occupancy[[i]],1,sum))
    p = apply(out$occupancy[[i]],2,sum) / N
    persist = numeric(R)
    persist[p != 0] = 1
    gamma_div[i] = sum(persist)
  }
  return(cbind(alpha_div, gamma_div))
}


