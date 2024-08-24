  // parameters for the initial state probability vector
  array[D] matrix[S - 1, K_i] beta_i_raw;  // first column is zero
  // parameters for the transition probabilities
  array[D, S] matrix[S - 1, K_s] beta_s_raw; // first row is zero
  // parameters for the emissions probabilities
  array[D] row_vector[S * (sum(M) - C) * K_o] beta_o_raw; // first row is zero
  matrix[D - 1, K_d] theta_raw;  // first column is zero
