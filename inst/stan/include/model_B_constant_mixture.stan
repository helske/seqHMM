  for(s in 1:S) {
    log_B[s, 1:M] = log_softmax(beta_o[d, s, , 1])';
    log_B[s, M + 1] = 0;
  }
