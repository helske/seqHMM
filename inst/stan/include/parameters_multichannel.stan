  // parameters for the initial state probability vector
  matrix[S - 1, K_i] beta_i_raw;  // first column is zero
  // parameters for the transition probabilities
  array[S] matrix[S - 1, K_s] beta_s_raw; // first row is zero
  // parameters for the emissions probabilities
  row_vector[S * (sum(n_M) - C) * K_o] beta_o_raw; // first row is zero
