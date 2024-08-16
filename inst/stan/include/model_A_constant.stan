    for(s in 1:S) {
      log_A[s, ] = log_softmax(beta_s[s, , 1])';
    }
