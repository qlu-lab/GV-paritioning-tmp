matrixV = function(GRM, v1_n, N) {
  # '''
  # generate (v1_n)th matrix V based on GRM 
  # Vs are symmetric matrices
  # GRM : 2N*2N
  # v1_n : index of wanted matrix
  # N : dimension of matrix
  # return: list that contains four blocks for (v1_n)th matrix V, each block is either 0, or a N-by-N matrix
  # '''
  V = list()
  V[[1]] = 0 ## top left block of V
  V[[2]] = 0 ## top right block
  V[[3]] = 0 ## bottom left block
  V[[4]] = 0 ## bottom right block
  if (v1_n == 1) {
    index = 1:N
    V[[1]] = GRM[index, index] - diag(1, N)
  } else if (v1_n == 2) {
    index = (N + 1):(2 * N)
    V[[1]] = 4 * GRM[index, index] - diag(1, N) * 2
  } else if (v1_n == 3) {
    M = matrix(c(0, 1, 1, 0), 2, 2)
    V[[1]] = kronecker(diag(1, N / 2), M)
  } else if (v1_n == 12) {
    index1 = 1:N
    index2 = (N + 1):(2 * N)
    V[[1]] = 2 * (GRM[index1, index2] + GRM[index2, index1]) - diag(1, N) * 2
  } else if (v1_n == 4) {
    index = 1:N
    V[[2]] = GRM[index, index]
    V[[3]] = V[[2]]
  } else if (v1_n == 5) {
    index = (N + 1):(2 * N)
    V[[2]] = GRM[index, index] * 4
    V[[3]] = V[[2]]
  } else if (v1_n == 6) {
    V[[2]] = diag(1, N)
    V[[3]] = V[[2]]
  } else if (v1_n == 7) {
    M = matrix(c(0, 1, 1, 0), 2, 2)
    V[[2]] = kronecker(diag(1, N / 2), M)
    V[[3]] = V[[2]]
  } else if (v1_n == 11) {
    index1 = 1:N
    index2 = (N + 1):(2 * N)
    V[[2]] = 2 * GRM[index1, index2]
    V[[3]] = 2 * GRM[index2, index1]
  } else if (v1_n == 14) {
    index1 = 1:N
    index2 = (N + 1):(2 * N)
    V[[3]] = 2 * GRM[index1, index2]
    V[[2]] = 2 * GRM[index2, index1]
  } else if (v1_n == 8) {
    index = 1:N
    V[[4]] = GRM[index, index] - diag(1, N)
  } else if (v1_n == 9) {
    index = (N + 1):(2 * N)
    V[[4]] = 4 * GRM[index, index] - diag(1, N) * 2
  } else if (v1_n == 10) {
    M = matrix(c(0, 1, 1, 0), 2, 2)
    V[[4]] = kronecker(diag(1, N / 2), M)
  } else if (v1_n == 13) {
    index1 = 1:N
    index2 = (N + 1):(2 * N)
    V[[4]] = 2 * (GRM[index1, index2] + GRM[index2, index1]) - diag(1, N) * 2
  }
  return(V)
}

