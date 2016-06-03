computeVP <- function(out){
  res = list()
  MEMout <- MEM_distMatrix(dist(out$XY))
  MEM.eigvec <- MEMout$MEM_IBD[,which(MEMout$eigenvaluesIBD>0)]
  
  speciesMatrix <- out$occupancy[[length(out$occupancy)]]
  res$bc <- RDA(speciesMatrix, MEM.eigvec)$Rsq_Adj
  
  
  res$ab = RDA(speciesMatrix, out$E)$Rsq_Adj
  
  res$abc = RDA(speciesMatrix, cbind(MEM.eigvec, out$E))$Rsq_Adj
  res$a = res$abc - res$bc
  res$c = res$abc - res$ab
  res$b = res$abc - res$a - res$c
  res$d = 1 - res$abc
}