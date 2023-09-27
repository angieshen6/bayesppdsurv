## usethis namespace: start
#' @importFrom Rcpp evalCpp
#' @useDynLib BayesPPDSurv,.registration = TRUE
## usethis namespace: end
NULL

#' Bayesian sample size determination using the power and normalized power prior for survival data
#'
#' The \pkg{BayesPPDSurv} (Bayesian Power Prior Design for Survival Data) package provides two categories of functions:
#' functions for Bayesian power/type I error calculation and functions for model fitting.
#' 
#' @details
#'
#' We assume the first column of the covariate matrix is the treatment indicator,
#' and the corresponding parameter is \eqn{\beta_1}.
#' The hypotheses are given by
#' \deqn{H_0: \beta_1 \ge \delta} and \deqn{H_1: \beta_1 < \delta.}
#' 
#' This implementation of the method does not assume any particular distribution for the sampling priors.
#' The user is allowed to specify a vector or matrix of samples for \eqn{\theta} (matrix if \eqn{\theta} is of dimension >1) from any distribution, and the algorithm samples with replacement
#' from the vector or matrix at each iteration of data simulation. In order to accurately approximate a joint distribution
#' for multiple parameters, the number of iterations should be large (e.g., 10,000).
#'
#' @references Chen, Ming-Hui, et al. "Bayesian design of noninferiority trials for medical devices using historical data." Biometrics 67.3 (2011): 1163-1170.
#' @docType package
#' @name BayesPPDSurv-package
NULL
#> NULL
