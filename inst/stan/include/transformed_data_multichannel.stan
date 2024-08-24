  row_vector[K_i] zeros_row_K_i = rep_row_vector(0, K_i);
  row_vector[K_s] zeros_row_K_s = rep_row_vector(0, K_s);
  row_vector[K_o] zeros_row_K_o = rep_row_vector(0, K_o);
  matrix[max_M, K_o] zeros_mat_M_K_o = rep_matrix(0, max_M, K_o);
  vector[S] zeros_S = rep_vector(0, S);
  row_vector[max_M + 1] zeros_M = rep_row_vector(0, max_M + 1);
  array[C + 1] int cumsum_M = append_array({0}, cumulative_sum(M));
