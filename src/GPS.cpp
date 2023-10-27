// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

void print(const char* str,
           const Eigen::MatrixXd& m) {
  Rprintf("%s\n", str);
  for (int i = 0; i < m.rows(); i ++) {
    for (int j = 0; j < m.cols(); j++ ) {
      Rprintf(" %g", m(i,j));
    }
    Rprintf("\n");
  }
}

void print(const char* str,
           const Eigen::VectorXd& m) {
  Rprintf("%s\n", str);
  for (int i = 0; i < m.size(); i ++) {
    Rprintf(" %g", m(i));
  }
  Rprintf("\n");
}

//' Rewrite lasso.sum.ess
//'
//' @param bVec genetic effect estimates
//' @param sVec standard error of genetic effect estimates;
//' @param r2Mat residual errors;
//' @param n sample size;
//' @param group 1-based group indicator
//' @param lambda l1 penalty parameter
//' @param alpha l2 penalty parameter
//' @param initVec initial value of beta
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
  // print("betaVec = ", betaVec);
  // beta0.vec <- rep(0,length(z.vec));
  const int vecLen = uVec.size();
  Eigen::VectorXd beta0Vec = Eigen::VectorXd::Zero(vecLen);
  
  double x_j_times_r, x_j_times_x;
  double old_beta_jj;
  double n_times_alpha_group_jj;
  
  // print("betaVec = ", betaVec);
  
  int iter  = 0;
  bool converged = false;
  // while(sum(abs(beta.vec-beta0.vec))>1e-5) {
  while ( !converged && iter < maxIter) {
    converged = (beta0Vec - betaVec).array().square().sum() <= 1e-5;
    // print("betaVec = ", betaVec);
    ++iter;
    // beta0.vec <- beta.vec;
    beta0Vec = betaVec;
    // for(jj in 1:length(z.vec)) {
    for (int jj = 0; jj < vecLen; ++jj) {
      // ##x_j * r_j = x_j (y - x_j' * beta_j' ) ;
      //x.j.times.r <- u.vec[jj]-sum(cov.mat[jj,-jj]*beta.vec[-jj])+2*n*eta*b.cc.vec[jj];
      old_beta_jj = betaVec(jj);
      betaVec(jj) = 0;
      x_j_times_r = uVec(jj) - (xTx.row(jj) * betaVec).sum() + 2.0 * n * eta(jj) * bccVec(jj);
      // Rprintf("x_j_times_r = %g\n", x_j_times_r);
      
      betaVec(jj) = old_beta_jj;
      
      // x.j.times.x <- cov.mat[jj,jj];
      x_j_times_x = xTx(jj, jj);
      // Rprintf("x_j_times_x = %g\n", x_j_times_x);
      // if(x.j.times.r > n*alpha[group[jj]] )
      //   beta.vec[jj] <- (x.j.times.r-n*(alpha[group[jj]]))/(x.j.times.x + 2*n*lambda[group[jj]]);
      // if(x.j.times.r < -n*alpha[group[jj]] )
      //   beta.vec[jj] <- (x.j.times.r+n*(alpha[group[jj]]))/(x.j.times.x + 2*n*lambda[group[jj]]);
      // if(x.j.times.r < n*alpha[group[jj]] & x.j.times.r > -n*alpha[group[jj]])
      //   beta.vec[jj] <- 0;
      n_times_alpha_group_jj = alpha( group(jj) - 1) * n;
      // Rprintf("n_times_alpha_group_jj = %g\n", n_times_alpha_group_jj);
      
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
  // return(beta.vec);
  return List::create(Named("beta") = betaVec,
                      Named("iteration") = iter,
                      Named("isConverged") = converged);
}





















// List fast_lasso_sum_ess(const Eigen::VectorXd& bVec,
//                         const Eigen::VectorXd& sVec,
//                         const Eigen::MatrixXd& r2Mat,
//                         int n,
//                         const Eigen::VectorXi& group,
//                         const Eigen::VectorXd& lambda,
//                         const Eigen::VectorXd& alpha,
//                         const Eigen::VectorXd& initVec,
//                         const int maxIter = 100) {
// 
//   // z.vec <- b.vec/s.vec;
//   Eigen::VectorXd zVec = bVec.array() / sVec.array();
//   // u.vec <- b.vec/s.vec^2;
//   Eigen::VectorXd uVec = bVec.array() / sVec.array().square();
//   // v.vec <- 1/s.vec;
//   Eigen::VectorXd vVec = sVec.array().inverse();
// 
//   // cov.mat <- cor2cov(r2.mat,v.vec);
//   Eigen::MatrixXd covMat = r2Mat;
//   covMat.array().colwise() *= vVec.array();
//   covMat.array().rowwise() *= vVec.transpose().array();
//   // Eigen::MatrixXd covMat = vVec.asDiagonal() * r2Mat * vVec.asDiagonal();
// 
//   // print("vVec = ", vVec);
//   // print("r2Mat = ", r2Mat);
//   // print("covMat = ", covMat);
//   // // #beta.vec <- ginv(cov.mat)%*%u.vec;
//   // beta.vec <- b.vec;
//   Eigen::VectorXd betaVec = initVec;
//   // print("betaVec = ", betaVec);
//   // beta0.vec <- rep(0,length(z.vec));
//   const int vecLen = bVec.size();
//   Eigen::VectorXd beta0Vec = Eigen::VectorXd::Zero(vecLen);
// 
//   double x_j_times_r, x_j_times_x;
//   double old_beta_jj;
//   double n_times_alpha_group_jj;
// 
//   // print("betaVec = ", betaVec);
// 
//   int iter  = 0;
//   bool converged = false;
//   // while(sum(abs(beta.vec-beta0.vec))>1e-5) {
//   while ( !converged && iter < maxIter) {
//     converged = (beta0Vec - betaVec).array().square().sum() <= 1e-5;
//     // print("betaVec = ", betaVec);
//     ++iter;
//     // beta0.vec <- beta.vec;
//     beta0Vec = betaVec;
//     // for(jj in 1:length(z.vec)) {
//     for (int jj = 0; jj < vecLen; ++jj) {
//       // ##x_j * r_j = x_j (y - x_j' * beta_j' ) ;
//       //x.j.times.r <- u.vec[jj]-sum(cov.mat[jj,-jj]*beta.vec[-jj]);
//       old_beta_jj = betaVec(jj);
//       betaVec(jj) = 0;
//       x_j_times_r = uVec(jj) - (covMat.row(jj) * betaVec).sum();
//       // Rprintf("x_j_times_r = %g\n", x_j_times_r);
// 
//       betaVec(jj) = old_beta_jj;
// 
//       // x.j.times.x <- cov.mat[jj,jj];
//       x_j_times_x = covMat(jj, jj);
//       // Rprintf("x_j_times_x = %g\n", x_j_times_x);
//       // if(x.j.times.r > n*alpha[group[jj]] )
//       //   beta.vec[jj] <- (x.j.times.r-n*(alpha[group[jj]]))/(x.j.times.x + 2*n*lambda[group[jj]]);
//       // if(x.j.times.r < -n*alpha[group[jj]] )
//       //   beta.vec[jj] <- (x.j.times.r+n*(alpha[group[jj]]))/(x.j.times.x + 2*n*lambda[group[jj]]);
//       // if(x.j.times.r < n*alpha[group[jj]] & x.j.times.r > -n*alpha[group[jj]])
//       //   beta.vec[jj] <- 0;
//       n_times_alpha_group_jj = alpha( group(jj) - 1) * n;
//       // Rprintf("n_times_alpha_group_jj = %g\n", n_times_alpha_group_jj);
// 
//       if (x_j_times_r > n_times_alpha_group_jj) {
//         betaVec(jj) = (x_j_times_r - n_times_alpha_group_jj) / (x_j_times_x + 2.0 * n * lambda( group(jj) - 1) );
//       } else if (x_j_times_r <  - n_times_alpha_group_jj) {
//         betaVec(jj) = (x_j_times_r + n_times_alpha_group_jj) / (x_j_times_x + 2.0 * n * lambda( group(jj) - 1) );
//       } else if (x_j_times_r < n_times_alpha_group_jj && x_j_times_r > -n_times_alpha_group_jj){
//         betaVec(jj) = 0;
//       } else {
//         REprintf("something wrong!\n");
//       }
//     }
//   }
//   if (iter == maxIter) {
//     REprintf("max iteration reached!!\n");
//   }
//   // return(beta.vec);
//   return List::create(Named("beta") = betaVec,
//                       Named("iteration") = iter,
//                       Named("isConverged") = converged);
// }

// sourceCpp('rcppeigen_hello_world.cpp')
