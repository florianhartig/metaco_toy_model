#' Plots an overview of the model output
#' @param out model output
plotOutput <- function(out, folder = "."){
  
  nSteps = dim(out$occupancy)[1]
  
  stats = computeOccupancyStatistics(out)
  
  par(mfrow = c(2,2))
  
  
  #pdf(file = paste(folder, "/spatial.pdf", sep = ""))
  
  plot(out$XY[,1], out$XY[,2], col=plotrix::color.scale(out$E[,1],c(0,1,1),c(1,1,0),0), main = "Environment (col) / occupancy \n (size) / turnover (ring)", pch=16,cex= 2* stats$occupancy$spatial / max(stats$occupancy$spatial), xlab = "Longitude", ylab = "Lattitude")
  points(out$XY[,1], out$XY[,2], cex= 2* stats$occupancy$spatial / max(stats$occupancy$spatial) , pch = 1)
  points(out$XY[,1], out$XY[,2], cex= 2* stats$occupancy$spatial / max(stats$occupancy$spatial) * stats$turnover$spatial / mean(stats$turnover$spatial) , pch = 1)
  
  #dev.off()

  
  
  barplot(log10(c(mean(stats$richness$alpha), mean(stats$richness$gamma), stats$occupancy$mean,stats$occupancy$meanAny, stats$turnover$mean)), names.arg = c("alpha", "gamma", "occupancy", "occupancyAny", "turnover"), las= 2, horiz  = T, main = "Occupacy statistics", xlab = "Value [Log10]")
  
  
  plot(stats$occupancy$species, stats$turnover$species / stats$occupancy$species, main = "Species Occupancy / Turnover / O")
  
  VP = computeVP(out)
  
  barplot(unlist(VP), names.arg = names(VP), las = 2, main = "Variation Partitioning")
  
  
#   #dev.new(width = 8, height = 6)
#   par(mar = c(5,6,2,1))
#   plot(c(1:nSteps), alpha_div, type = "l", xlab = "Time", ylab = "Diversity", cex.lab = 1.25, cex.axis = 1.25)
#   lines(loess.smooth(c(1:nSteps), alpha_div), col = "red")
#   title(main = "Local diversity")
#   
#   #dev.new(width = 8, height = 6)
#   par(mar = c(5,6,2,1))
#   plot(c(1:nSteps), gamma_div, type = "l", xlab = "Time", ylab = "Diversity", cex.lab = 1.25, cex.axis = 1.25)
#   lines(loess.smooth(c(1:nSteps), gamma_div), col = "red")
#   
#   title(main = "Regional diversity")
  
  
  #plot(out$colonization, type = "l", main = "Colonization")
  #lines(loess.smooth(c(1:nSteps), out$colonization), col = "red")
  
  return(list(occupancyStatistics = stats, VP = VP))
  
  
}

#' Computes occupancy statistics 
#' @param out result of a simulation
#' @return a list with richness (alpha, gamma), occupancy (temporal, spatial, species) and turnover (temporal, spatial, species) 
computeOccupancyStatistics <- function(out, start = 1){
  # Compute the number of species per time step
  nSteps = dim(out$occupancy)[1] - start + 1
  nPatches = dim(out$occupancy)[2]
  nSpecies = dim(out$occupancy)[3]
  
  timesteps = (start:(start + nSteps - 1))
  
  alpha_div = numeric(nSteps)
  gamma_div = numeric(nSteps)
  occupancy = list(temporal = numeric(nSteps), temporalAny = numeric(nSteps), spatial = numeric(nPatches), spatialAny = numeric(nPatches), species = numeric(nSpecies))
  turnover = list(temporal = numeric(nSteps), temporalAny = numeric(nSteps),spatial = numeric(nPatches), spatialAny = numeric(nPatches),species = numeric(nSpecies))
  
  
  for(i in start:(start + nSteps - 1)) {
    
    # richness
    
    alpha_div[i] = mean(apply(out$occupancy[i,,],1,sum))
    p = apply(out$occupancy[i,,],2,sum) / N
    persist = numeric(R)
    persist[p != 0] = 1
    gamma_div[i] = sum(persist)
    
    # turnover and occupancy
    occupancy$temporal[i] = mean(out$occupancy[i,,])
    occupancy$spatial = occupancy$spatial + apply(out$occupancy[i,,],1,mean)
    anyOccupancyPatch = apply(out$occupancy[i,,],1,max)
    occupancy$spatialAny = occupancy$spatialAny + anyOccupancyPatch
    occupancy$temporalAny[i] = mean(anyOccupancyPatch)
    
    occupancy$species = occupancy$species + apply(out$occupancy[i,,], 2,mean)
    
    if(i>1){
      turn = abs(out$occupancy[i,,] - out$occupancy[i-1,,])
      
      turnover$temporal[i] =  mean(turn)
      
      turnover$spatial = turnover$spatial + apply(turn,1,mean)
      anyturnoverPatch = apply(turn,1,max)
      turnover$spatialAny = turnover$spatialAny + anyturnoverPatch
      turnover$temporalAny[i] = mean(anyturnoverPatch)
      
      turnover$species = turnover$species + apply(turn, 2,mean)
      
    }
  }
  
  occupancy$spatial = occupancy$spatial / nSteps
  occupancy$spatialAny = occupancy$spatialAny / nSteps
  occupancy$species = occupancy$species / nSteps
  
  turnover$spatial = turnover$spatial / nSteps
  turnover$spatialAny = turnover$spatialAny / nSteps
  turnover$species = turnover$species / nSteps
  
  occupancy$mean = mean(occupancy$temporal)
  occupancy$meanAny = mean(occupancy$temporalAny)
  
  turnover$mean = mean(turnover$temporal)
  turnover$meanAny = mean(turnover$temporalAny)
  
  #$temporal[round(nSteps/2):round(nSteps/2)]
  
  return(list( richness = list(alpha = alpha_div, gamma = gamma_div), occupancy =  occupancy, turnover = turnover))
}

# stats = computeOccupancyStatistics(out)
# 
# par(mfrow = c(3,3))
# 
# 
# plot(stats$occupancy$temporal)
# plot(stats$occupancy$temporalAny)
# barplot(stats$occupancy$spatial)
# barplot(stats$occupancy$spatialAny)
# barplot(stats$occupancy$species)
# 
# plot(stats$turnover$temporal)
# plot(stats$turnover$temporalAny)
# barplot(stats$turnover$spatial)
# barplot(stats$turnover$spatialAny)
# barplot(stats$turnover$species)



