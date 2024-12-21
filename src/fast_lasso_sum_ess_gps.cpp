// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
//' @title core function for fast coordinate descent of GPS loss function
//' @param uVec u statistics (x'y)
//' @param xTx x_transpose_x matrix calculated from reference panel
//' @param n sample size;
//' @param group 1-based group indicator that indicates which alpha and lambda to use
//' @param alpha l1 penalty parameter options, match unique group id
//' @param lambda l2 penalty parameter options, match unique group id
//' @param eta l2 penalty parameter for deviation from prior (correlated) effect size estimates
//' @param bccVec weight of PRS for a correlated trait to serve as a prior
//' @param initVec initial value of beta (joint effect)
//' @param maxIter maximum number of iterations
//' @return a list where `beta` is the fitted regression coefficient vector, and `iteration` is the actual iteration.
// [[Rcpp::export]]
List fast_lasso_sum_ess_gps(const Eigen::VectorXd& uVec,
                            const Eigen::MatrixXd& xTx,
                            int n,
                            const Eigen::VectorXi& group,
                            const Eigen::VectorXd& lambda,
                            const Eigen::VectorXd& alpha,
                            const Eigen::VectorXd& eta,
                            const Eigen::VectorXd& bccVec,
                            const Eigen::VectorXd& initVec,
                            const int maxIter = 100) {
  // Initialize the joint effect estimate
  Eigen::VectorXd betaVec = initVec;
  const int vecLen = uVec.size();
  Eigen::VectorXd beta0Vec = Eigen::VectorXd::Zero(vecLen);

  double x_j_times_r, x_j_times_x;
  double old_beta_jj;
  double n_times_alpha_group_jj;
  int iter  = 0;
  bool converged = false;
  while ( !converged && iter < maxIter) {
    converged = (beta0Vec - betaVec).array().square().sum() <= 1e-5;
    ++iter;
    beta0Vec = betaVec;
    for (int jj = 0; jj < vecLen; ++jj) {
      old_beta_jj = betaVec(jj);
      betaVec(jj) = 0;
      x_j_times_r = uVec(jj) - (xTx.row(jj) * betaVec).sum() + 2.0 * n * eta(jj) * bccVec(jj);
      betaVec(jj) = old_beta_jj;
      x_j_times_x = xTx(jj, jj);
      n_times_alpha_group_jj = alpha( group(jj) - 1) * n;
      if (x_j_times_r > n_times_alpha_group_jj) {
        betaVec(jj) = (x_j_times_r - n_times_alpha_group_jj) / (x_j_times_x + 2.0 * n * lambda( group(jj) - 1) + 2.0 * n * eta(jj));
      } else if (x_j_times_r <  - n_times_alpha_group_jj) {
        betaVec(jj) = (x_j_times_r + n_times_alpha_group_jj) / (x_j_times_x + 2.0 * n * lambda( group(jj) - 1) + 2.0 * n * eta(jj));
      } else if (x_j_times_r < n_times_alpha_group_jj && x_j_times_r > -n_times_alpha_group_jj){
        betaVec(jj) = 0;
      } else {
        REprintf("something wrong!\n");
      }
    }
  }
  if (iter == maxIter) {
    REprintf("max iteration reached!!\n");
  }
  return List::create(Named("beta") = betaVec,
                      Named("iteration") = iter,
                      Named("isConverged") = converged);
}
