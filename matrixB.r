

matrixB = function(GRM, y1, y2, n_para, N, n_cores) {
  # '''
  # matrix on the right hand side of linear system
  # matrix B = yVy-trace(V)
  # GRM : 2N-by-2N GRM matrix 
  # n_para : number of parameters to estimate
  # N : dimension of matrix
  # n_cores: cores used for processing
  # '''
  A = matrix(0, n_para, 1)
  rowlist = c(1:n_para)
  B = c()
  for (row in 1:n_para) {
    print(row)
    
    V = matrixV(GRM, row, N)
    term1 = 0
    term2 = 0
    for (i in 1:4) {
      if (is.null(dim(V[[i]]))) {
        next
      }
      if (i == 1) {
        term1 = term1 + t(y1) %*% V[[1]] %*% y1
      } else if (i == 2) {
        term1 = term1 + t(y1) %*% V[[2]] %*% y2
      } else if (i == 3) {
        term1 = term1 + t(y2) %*% V[[3]] %*% y1
      } else if (i == 4) {
        term1 = term1 + t(y2) %*% V[[4]] %*% y2
      }
      term2 = term2 + sum(diag(V[[i]]))
    }
    B = c(B, term1 - term2)
    
  }
  
  return(B)
}
