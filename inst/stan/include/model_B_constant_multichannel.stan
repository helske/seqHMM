    for (c in 1:C) {
      for (s in 1:S) {
        log_B[c, s, 1:M[c]] = log_softmax(beta_o[c, s, 1:M[c], 1])';
        log_B[c, s, M[c] + 1] = 0;
      }
    }
