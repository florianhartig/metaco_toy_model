computeVP <- function(out){
  
  nSteps = dim(out$occupancy)[1]
  
  speciesMatrix <- out$occupancy[nSteps,,]
  
  
  absentSpecies = apply(speciesMatrix, 2, sum) == 0
  emptySites = apply(speciesMatrix, 1, sum) == 0
  
  absentSpecies = rep(F, D)
  emptySites = rep(F, R)
  
  speciesMatrix = speciesMatrix[,!absentSpecies]
  speciesMatrix = speciesMatrix[!emptySites, ]
  
  res = list()
  MEMout <- MEM_distMatrix(dist(out$XY[!emptySites, ]))
  MEM.eigvec <- MEMout$MEM_IBD[,which(MEMout$eigenvaluesIBD>0)]
  
  res$bc <- RDA(speciesMatrix, MEM.eigvec)$Rsq_Adj
  
  
  res$ab = RDA(speciesMatrix, out$E[!emptySites,])$Rsq_Adj
  
  res$abc = RDA(speciesMatrix, cbind(MEM.eigvec, out$E[!emptySites,]))$Rsq_Adj
  res$a = res$abc - res$bc
  res$c = res$abc - res$ab
  res$b = res$abc - res$a - res$c
  res$d = 1 - res$abc
  
  return(res)
}