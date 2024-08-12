array[C, S] matrix[M, K_o] beta_o;
{
  array[K_o] int k_range = linspaced_int_array(K_o, 1, K_o);
  row_vector[K_o] zeros = rep_row_vector(0, K_o);
  for (c in 1:C) {
    for (s in 1:S) {
      beta_o[c, s, 1, ] = zeros;
      for (m in 2:M) {
        array[K_o] int idx = add(k_range, (c - 1) * S * (M - 1) * K_o + (s - 1) * (M - 1) * K_o + (m - 2) * K_o);
        beta_o[c, s, m, ] = beta_o_raw[idx];
      }
    }
  }
}
