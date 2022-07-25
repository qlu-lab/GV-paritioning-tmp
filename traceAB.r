
traceAB = function(V1, V2, n_cores = 1) {
  # '''
  # A and B are both blockwise matrices with four blocks like [[A1,A2],[A3,A4]]
  # trace(AB) = trace(A1B1+A2B3+A3B2+A4B4)
  # V1,V2 : output from matrixV()
  # '''
  
  i_v <- c(1, 2, 3, 4)
  j_v <- c(1, 3, 2, 4)
  
  result = c()
  for (k in 1:4) {
    i = i_v[k]
    j = j_v[k]
    if ((!is.null(dim(V1[[i]]))) & (!is.null(dim(V2[[j]])))) {
      result = c(result, sum(V1[[i]] * V2[[j]]))
    } else{
      result = c(result, 0)
    }
  }
  
  tr = sum(result)
  return(tr)
}
