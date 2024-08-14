matrix[S, K_i] beta_i = append_row(zeros_row_K_i, beta_i_raw);
array[S] matrix[S, K_s] beta_s;
array[C, S] matrix[M, K_o] beta_o;
for(s in 1:S) {
  beta_s[s] = append_row(zeros_row_K_s, beta_s_raw[s]);
}
{
  array[K_o] int k_range = linspaced_int_array(K_o, 1, K_o);
  for (c in 1:C) {
    for (s in 1:S) {
      beta_o[c, s, 1, ] = zeros_row_K_o;
      for (m in 2:M) {
        array[K_o] int idx = add(k_range, (c - 1) * S * (M - 1) * K_o + (s - 1) * (M - 1) * K_o + (m - 2) * K_o);
        beta_o[c, s, m, ] = beta_o_raw[idx];
      }
    }
  }
}
