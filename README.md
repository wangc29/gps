Chen Wang
2023-12-30

GPS
=====

**G**enetic **P**rediction **S**core

GPS is a penalized regression based method that can integrate large scale case-control study and electronic health record (EHR) based biobanks to construct polygenic risk score models for disease progression risk from preclinical stage with improved accuracy.

GPS is tested on R 3.6+ for Linux and macOS systems and depends on the following packages: Rcpp,RcppEigen

GPS can be installed from our github repo by utilizing the devtools library. Installation only takes a few minutes.

Installation
------------

    library(devtools)
    install_github("wangc29/gps")

Fit GPS models
--------------

```{r}
eta_vec<-c(0.5,1.25,2,2.75,3.5)
NN<-10000
gps_res<-gps_grid(lambda_max = max(abs(beta_est)),min_ratio = 0.01,nlambda = 30,
                  eta_vec = eta_vec,beta_est = beta_est,se_est = beta_se,
                  ld_cor_list = ld_cor_list,sample_size = NN,
                  weight_prior = beta_cc_est,init_vec = NULL,max_iter = 1000)
```
where

- eta_vec: tuning parameters for the contribution of prior estimate of effect size from case-control study
- lambda: maximum value of tuning parameter for shrinkage
- beta_est: marginal association effect size estimate vector for progression trait
- se_est: marginal standard error estumates of beta_est
- ld_cor_list: a list of LD correlation matrix.
- sample_size: sample size of progression trait cohort
- weight_prior: joint effect estimates (used as prior) from case-control studies  for corresponding variants in beta_est.

The function returns a matrix containing joint effect size estimates for progression trait. Each column respresents joint effect size estimate vector that corresponds to a specific configuration of tuning parameters. The matrix can be retrieved by,

```R
beta_joint_est_matrix = gps_res[["beta_grid"]]
```

The corresponding tuning parameters can be retrieved by

```R
beta_joint_est_matrix = gps_res[["param_grid"]]
```



