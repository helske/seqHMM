// Get correct index for W_A and W_B based on response values
#include "get_W_idx.h"

arma::uword get_W_idx(const arma::uvec& y, const arma::uvec& M) {
  arma::uword index = 0;
  arma::uword multiplier = 1;
  arma::uword C = y.n_elem;
  
  for (arma::uword c = 0; c < C; ++c) {
    index += y[c] * multiplier;
    multiplier *= M[c];
  }
  return index;
}
