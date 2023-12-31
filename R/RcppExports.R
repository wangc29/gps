# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Rewrite lasso.sum.ess
#'
#' @param bVec genetic effect estimates
#' @param sVec standard error of genetic effect estimates;
#' @param r2Mat residual errors;
#' @param n sample size;
#' @param group 1-based group indicator
#' @param lambda l1 penalty parameter
#' @param alpha l2 penalty parameter
#' @param initVec initial value of beta
#' @param maxIter maximum number of iterations
#' @return a list where `beta` is the fitted regression coefficient vector, and `iteration` is the actual iteration.
fast_lasso_sum_ess_gps <- function(uVec, xTx, n, group, lambda, alpha, eta, bccVec, initVec, maxIter = 100L) {
    .Call(`_GPS_fast_lasso_sum_ess_gps`, uVec, xTx, n, group, lambda, alpha, eta, bccVec, initVec, maxIter)
}

