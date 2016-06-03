MEM_distMatrix<-function(Dist_matrix,matrix_type="IBD") 
{
  n <- nrow(Dist_matrix)
  
  if (matrix_type=="IBD") {
    ####################################
    # IBD treatment (Euclidean matrix) #
    ####################################
    IBD_matrix <- as.matrix(Dist_matrix)
    # truncate matrix
    require(vegan)
    spanning <- vegan::spantree(IBD_matrix)
    threshh <- max(spanning$dist)
    IBD_matrix_truncated <- IBD_matrix
    IBD_matrix_truncated[IBD_matrix_truncated > threshh] <- 4 * threshh
    # note that IBD_matrix_truncated has a main diagonal of zero so that is congrouent with the MEM framework as in Dray et al. (2006)
    # double center truncated matrix
    diag(IBD_matrix_truncated) <- 4 * threshh
    row.wt = rep(1, nrow(IBD_matrix_truncated))
    col.wt = rep(1, ncol(IBD_matrix_truncated))
    st <- sum(col.wt)
    sr <- sum(row.wt)
    row.wt <- row.wt/sr
    col.wt <- col.wt/st
    Centered_Matrix <- -0.5*(IBD_matrix_truncated*IBD_matrix_truncated)
    row.mean <- apply(row.wt * Centered_Matrix, 2, sum)
    col.mean <- apply(col.wt *t(Centered_Matrix), 2, sum)
    col.mean <- col.mean - sum(row.mean * col.wt)
    Centered_Matrix <- sweep(Centered_Matrix, 2, row.mean)
    Centered_Matrix <- t(sweep(t(Centered_Matrix), 2, col.mean)) 
    D.eig <- eigen(Centered_Matrix, symmetric = TRUE)
    epsilon <- sqrt(.Machine$double.eps)
    Zero_eigenvalue<-which(abs(D.eig$values) < epsilon) # due to centering, one eigenvector has zero variance
    eigenvalues_IBD <- D.eig$values[-Zero_eigenvalue]
    eigenvectors_IBD <- subset(D.eig$vectors, select= -Zero_eigenvalue)
    weight <- sqrt(abs(eigenvalues_IBD))
    eigenvectors_IBD <- eigenvectors_IBD%*% diag(weight)
    result <- list (eigenvaluesIBD = eigenvalues_IBD, MEM_IBD=eigenvectors_IBD)
  }  
  
  if (matrix_type=="IBR") {
    ########################################
    # IBR treatment (non-Euclidean matrix) #
    ########################################
    IBR_matrix <- as.matrix(Dist_matrix)
    # the IBR matrix is already non-Euclidean by design and it won't be truncated
    row.wt = rep(1, nrow(IBR_matrix))
    col.wt = rep(1, ncol(IBR_matrix))
    st <- sum(col.wt)
    sr <- sum(row.wt)
    row.wt <- row.wt/sr
    col.wt <- col.wt/st
    Centered_Matrix <- -0.5*(IBR_matrix*IBR_matrix)
    row.mean <- apply(row.wt * Centered_Matrix, 2, sum)
    col.mean <- apply(col.wt *t(Centered_Matrix), 2, sum)
    col.mean <- col.mean - sum(row.mean * col.wt)
    Centered_Matrix <- sweep(Centered_Matrix, 2, row.mean)
    Centered_Matrix <- t(sweep(t(Centered_Matrix), 2, col.mean)) 
    D.eig <- eigen(Centered_Matrix, symmetric = TRUE)
    epsilon <- sqrt(.Machine$double.eps)
    Zero_eigenvalue<-which(abs(D.eig$values) < epsilon) # due to centering, one eigenvector has zero variance
    eigenvalues_IBR <- D.eig$values[-Zero_eigenvalue]
    eigenvectors_IBR <- subset(D.eig$vectors, select= -Zero_eigenvalue)
    weight <- sqrt(abs(eigenvalues_IBR))
    eigenvectors_IBR <- eigenvectors_IBR%*% diag(weight)
    result <- list (eigenvaluesIBR = eigenvalues_IBR, MEM_IBR=eigenvectors_IBR)
  }
  
  return(result)
}




