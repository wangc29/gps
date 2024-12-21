#' cor2cov function
#' @param ld.cor.mat a LD correlation matrix
#' @param v.vec a vector of v statistic
#' @return the approximated X'X from summary statistics and LD.
cor2cov <- function(ld.cor.mat, v.vec) {
  if (length(v.vec)==1) {
    v.vec<-rep(v.vec,ncol(ld.cor.mat))
  }
  diag(v.vec) %*% ld.cor.mat %*% diag(v.vec)
}

#' fit a single LD block
#'
#' @return a vector
gps_fit_single<-function(beta_est,se_est,ld_cor_mt,sample_size,anno_group,weight_prior,lambda,pi,phi,eta,init_vec=NULL,max_iter=1000,snp_id=NULL){
  ##initial joint effect est
  if (is.null(init_vec)) {
    init_vec<-beta_est
  }
  lambda_option<-c(lambda,lambda*phi)
  ##Now considering the mixing param of L1 and L2 penalty
  lambda_L1<-lambda_option*pi
  lambda_L2<-lambda_option*(1-pi)
  #weight prior: non-existing prior PRS weights are given as NA.
  ##If eta is given as 0, then shut down all pentalty to the prior
  if (eta>0) {
    eta_vec<-rep(eta,length(weight_prior))
    eta_vec[which(is.na(weight_prior))]<-0
    weight_prior[which(is.na(weight_prior))]<-0
  } else {
    eta_vec<-rep(0,length(weight_prior))
    weight_prior<-rep(0,length(weight_prior))
  }
  ##Calculate x'y and x'x
  ustat<- beta_est/(se_est)^2 ## x'y
  vstat<- 1/se_est
  x_transpose_x<-cor2cov(ld.cor.mat = ld_cor_mt,v.vec = vstat)
  ##Run core code of fitting GPS model
  res_lasso_gps<-fast_lasso_sum_ess_gps(uVec = ustat,xTx = x_transpose_x,n = sample_size,group = anno_group,
                                        lambda = lambda_L2,alpha = lambda_L1,eta = eta_vec,
                                        bccVec = weight_prior,initVec = init_vec ,maxIter = max_iter)
  res_lasso_gps[["snp_id"]] <- snp_id

  return(res_lasso_gps)
}

#' fit across LD blocks
#'
gps_fit_global<-function(beta_est,se_est,ld_cor_list,sample_size,anno_group,weight_prior,lambda,pi,phi,eta,init_vec=NULL,max_iter,snp_id){
  nBlock<-length(ld_cor_list)
  # Fit by blocks
  block_end_ix<-c(0,cumsum(sapply(ld_cor_list, function(x){dim(x)[1]})))
  beta_by_blk<-list()
  snp_id_by_blk<-list()
  beta_single_by_blk<-list()
  se_single_by_blk<-list()
  iteration_by_blk<-rep(NA,nBlock)
  isConverged_by_blk<-rep(NA,nBlock)
  for (blk_num in 1:nBlock) {
    #cat("Running block: ",blk_num,"\n")
    blk_start<-block_end_ix[blk_num]+1
    blk_end<-block_end_ix[blk_num+1]
    beta_est_kk<-beta_est[blk_start:blk_end]
    beta_single_by_blk[[blk_num]]<-beta_est_kk
    se_est_kk<-se_est[blk_start:blk_end]
    se_single_by_blk[[blk_num]]<-se_est_kk
    ld_cor_mt_kk<-ld_cor_list[[blk_num]]
    anno_group_kk<-anno_group[blk_start:blk_end]
    weight_prior_kk<-weight_prior[blk_start:blk_end]
    init_vec_kk<-NULL
    if (!is.null(init_vec)){init_vec_kk<-init_vec[blk_start:blk_end]};
    snp_id_kk<-snp_id[blk_start,blk_end]
    res_gps_kk<-gps_fit_single(beta_est = beta_est_kk,se_est = se_est_kk,ld_cor_mt = ld_cor_mt_kk,sample_size = sample_size,
                               anno_group = anno_group_kk, weight_prior = weight_prior_kk,lambda = lambda ,pi = pi,phi = phi,eta = eta,
                               init_vec = init_vec_kk,max_iter = max_iter ,snp_id = snp_id_kk)
    beta_by_blk[[blk_num]]<-res_gps_kk[["beta"]]
    snp_id_by_blk[[blk_num]]<-res_gps_kk[["snp_id"]]
    iteration_by_blk[blk_num]<-res_gps_kk[["iteration"]]
    isConverged_by_blk[blk_num]<-res_gps_kk["isConverged"]
    #cat("Finished block: ",blk_num,"\n")
  }
  return(list(beta_by_blk=beta_by_blk,iteration_by_blk=iteration_by_blk,
              isConverged_by_blk=isConverged_by_blk,snp_id_by_blk=snp_id_by_blk,
              beta_single_by_blk=beta_single_by_blk,se_single_by_blk=se_single_by_blk))
}

#' log sequence of hyperparams
#'
seq.log <- function(from, to, length.out) {
  return(exp(seq(log(from),log(to),length.out = length.out)))
}

#' fit grid
#'
#' @param min_ratio the ratio between the minimum lambda (hyperparameter for shrinkage) and the maximum lambda
#' @param nlambda number of lambda values that are evenly spaced beween min of lambda and max of lambda
#' @param l1_ratio the contribution of l1 penalty for shrinkage.
#' @param eta_vec a vector of eta values to be used. eta controls the contribution of the prior PRS model weights.
#' @param beta_est a vector of marginal effect size estimates (assuming standardized genotype and phenotype)
#' @param se_est a vector of standard error estimates of marginal effect size
#' @param ld_cor_list a list of LD correlation matrix, each element is a LD correlation matrix of a LD block
#' @param sample_size a scalar of sample size
#' @param weight_prior a vector of the PRS weights of a correlated traits (e.g. PRS of the case-control phenotype).
#' The order of the elements should be the same as beta_est and se_est. If not available, please fill with NAs.
#' @param init_vec initial values used for coordinate descent algorithm.
#' @param max_iter num of maximum iterations used in coordinate descent
#' @return List of two elements: 1. a grid of joint effect size estimates. Each column corresponds to a specific combination of hyperparameters
#' 2. a dataframe containing corresponding combinations of hyperparameters
#' @export
gps_grid<-function(min_ratio=0.01,nlambda=30,l1_ratio=c(0,0.25,0.5,0.75,1),eta_vec,
                   beta_est,se_est,ld_cor_list,sample_size,weight_prior,
                   init_vec=NULL,max_iter=1000){
  u_stat<-beta_est/(se_est)^2
  PP<-length(u_stat)
  lambda_max_vec<-max(abs(u_stat/sample_size))/l1_ratio

  param_grid_list<-list()
  for(ii in 1:length(lambda_max_vec)){
    lambda_max_ii<-lambda_max_vec[ii]
    lambda_vec_ii<-seq.log(lambda_max_ii*min_ratio,lambda_max_ii,length.out = nlambda)
    param_grid_ii<-expand.grid(lambda_param=lambda_vec_ii,eta_param=eta_vec)
    param_grid_ii[,"l1_ratio"]<-l1_ratio[ii]
    param_grid_list[[ii]]<-param_grid_ii
  }
  param_grid <-do.call('rbind',param_grid_list)
  ord<-nrow(param_grid)
  ##Fit GPS model for each tuning param combination
  beta_grid<-list()
  for (comb_num in 1:ord) {
    lambda_curr<-param_grid$lambda_param[comb_num]
    eta_curr<-param_grid$eta_param[comb_num]
    l1_ratio_curr<-param_grid$l1_ratio[comb_num]
    ##Fit gps on training data
    gps_res<-gps_fit_global(beta_est = beta_est,se_est =se_est,
                            ld_cor_list = ld_cor_list, sample_size = sample_size,
                            anno_group = rep(1,PP),weight_prior = weight_prior,
                            lambda = lambda_curr,pi = l1_ratio_curr,phi = c(1),eta = eta_curr,init_vec = init_vec,
                            max_iter = max_iter,snp_id = NULL)
    beta_joint<-unlist(gps_res$beta_by_blk)
    beta_grid[[comb_num]]<-beta_joint
  }
  return(list("beta_grid"=do.call(cbind,beta_grid),"param_grid"=param_grid))
}
