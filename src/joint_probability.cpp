#include "joint_probability.h"
// [[Rcpp::export()]]
arma::vec joint_probability(const arma::field<arma::vec>& prob) {
  arma::vec joint = prob(prob.n_elem - 1);
  for (arma::sword i = prob.n_elem - 2; i >= 0; --i) {
    joint = arma::kron(joint, prob(i));
  }
  return joint;
}