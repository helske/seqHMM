functions {
  array[] compute_initial(vector alpha_raw, matrix beta_raw, array[] vector X) {
    int S = size(alpha_raw) + 1;
    int N = dims(X)[1];
    int K = dims(X)[2];
    vector[S] alpha = reverse(append_row(0, alpha_raw));
    matrix[S, K] beta = append_row(rep_row_vector(0, K), beta_raw);
    array[N] vector[S] Pi;
    
    for(i in 1:N) {
      Pi[i] = softmax(alpha + beta * X[i]);
    }
    return Pi;
    
  }
  array[,] matrix compute_transitions(array[] matrix beta_raw, array[] matrix X) {
    int S = size(beta_raw);
    int T = size(X);
    int K = dims(X)[2];
    int N = dims(X)[3];
    array[S] matrix[S, K] beta;
    array[N, T] matrix[S, S] A;
    for(s in 1:S) {
      beta[s] = append_row(rep_row_vector(0, K), beta_raw[s]);
    }
    for(i in 1:N) {
      for(t in 1:T) {
        for(s in 1:S) {
          A[i, t, s, ] = softmax(beta[s] * X[t, , i])';
        }
      }
    }
    return A;
  }
  
  array[,] matrix compute_emissions(array[] matrix beta_raw, array[] matrix X) {
    int S = size(beta_raw);
    int T = size(X);
    int K = dims(X)[2];
    int N = dims(X)[3];
    int M = dims(beta_raw)[2] + 1;
    array[S] matrix[M, K] beta;
    for(s in 1:S) {
      beta[s] = append_row(rep_row_vector(0, K), beta_raw[s]);
    }
    array[N, T] matrix[S, M] B;
    for(i in 1:N) {
      for(t in 1:T) {
        for(s in 1:S) {
          B[i, t, s, ] = softmax(beta[s] * X[t, , i])';
        }
      }
    }
    return B;
  }
}

