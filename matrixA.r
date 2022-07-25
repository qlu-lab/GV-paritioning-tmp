
matrixA = function(n_para, GRM, N, n_cores) {
  # '''
  # matrix on the left hand side of linear system to solve
  # matrixA is symmetric matrix, so just need to calculate lower triangle
  # GRM : 2N-by-2N GRM matrix 
  # n_para : number of parameters to estimate
  # '''
  v1 = c(1,
         1,
         1,
         1,
         2,
         2,
         2,
         3,
         3,
         12,
         4,
         4,
         4,
         4,
         4,
         4,
         5,
         5,
         5,
         5,
         5,
         6,
         6,
         6,
         6,
         7,
         7,
         7,
         11,
         11,
         14)
  v2 = c(1,
         2,
         3,
         12,
         2,
         3,
         12,
         3,
         12,
         12,
         4,
         5,
         6,
         7,
         11,
         14,
         5,
         6,
         7,
         11,
         14,
         6,
         7,
         11,
         14,
         7,
         11,
         14,
         11,
         14,
         14)
  A = matrix(0, n_para, n_para)
  for (row in 1:length(v1)) {
    # a = Sys.time()
    V1 = matrixV(GRM, v1[row], N)
    V2 = matrixV(GRM, v2[row], N)
    A[v1[row], v2[row]] = traceAB(V1, V2, n_cores)
    # print(Sys.time() - a)
  }
  A[8, 8] = A[1, 1]
  A[9, 9] = A[2, 2]
  A[10, 10] = A[3, 3]
  A[8, 9] = A[1, 2]
  A[8, 10] = A[1, 3]
  A[9, 10] = A[2, 3]
  A[8, 13] = A[1, 12]
  A[9, 13] = A[2, 12]
  A[10, 13] = A[3, 12]
  A[13, 13] = A[12, 12]
  A = A + t(A)
  diag(A) = diag(A) / 2
  return(A)
}