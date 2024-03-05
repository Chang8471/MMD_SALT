
#' Count needed to call signal
#'
#' @description Count needed to reach desired postL for each sample and gene combination
#'
#' @param postL_cutoff cutoff used to binarize postL
#' @param useparEst_itr_out boolean indicator to use output object from ``parEst_itr`` function
#' @param parEst_itr_out if ``useparEst_itr_out`` is T, then provide the output object from ``parEst_itr`` function
#' @param gridResolution number of grid nodes to generate for the step function for posterior likelihood
#' @param alphai_beta_g mandatory parameter if ``useparEst_itr_out`` is F, will be ignored otherwise. Matrix of signal mean, after log normal approximation. Same dimensions as ``Y_mtx``.
#' @param sigma2_1g mandatory parameter if ``useparEst_itr_out`` is F, will be ignored otherwise. Matrix of signal variance, after log normal approximation. Same dimensions as ``Y_mtx``.
#' @param mu_0i mandatory parameter if ``useparEst_itr_out`` is F, will be ignored otherwise. Matrix of background mean. Same dimensions as ``Y_mtx``.
#' @param sigma2_0i mandatory parameter if ``useparEst_itr_out`` is F, will be ignored otherwise. Matrix of background variance. Same dimensions as ``Y_mtx``.
#' @param pi_prior mandatory parameter if ``useparEst_itr_out`` is F, will be ignored otherwise. Matrix of prior likelihood of signal. Same dimensions as ``Y_mtx``.
#'
#' @return Matrix of posterior likelihood, same dimensions as ``alphai_beta_g``.
#' @export
#'
#' @examples
calc_count_cutoff = function(postL_cutoff, useparEst_itr_out,  parEst_itr_out=NULL,gridResolution=500,
                     alphai_beta_g=NULL,sigma2_1g = NULL,mu_0i=NULL,sigma2_0i=NULL,pi_prior=NULL){

  if(useparEst_itr_out){
    if(!is.null(parEst_itr_out$S2_hat)){
      alphai_beta_g = parEst_itr_out$m_hat
      sigma2_1g = parEst_itr_out$S2_hat
    }else{
      alphai_beta_g = parEst_itr_out$alphai_beta_g
      sigma2_1g = parEst_itr_out$sigma2_1g
    }
    mu_0i = parEst_itr_out$mu_0i
    sigma2_0i = parEst_itr_out$sigma2_0i
    pi_prior = parEst_itr_out$pi_prior
  }

  if(any(is.null(alphai_beta_g),is.null(sigma2_1g),is.null(mu_0i),is.null(sigma2_0i),is.null(pi_prior))) stop("argument missing")


  ####### correction
  count_cufoffs = sapply(1:(nrow(alphai_beta_g)*ncol(alphai_beta_g)), function(i){
    P_signal_func = function(x) dlnorm(x,alphai_beta_g[i],sqrt(sigma2_1g[i]))
    P_background_func = function(x) dnorm(x, mu_0i[i], sqrt(sigma2_0i[i]))
    pi_prior_tmp = pi_prior[i]

    postL_func = function(x, postL_cutoff){
      P_signal_func(x)*pi_prior_tmp/(P_signal_func(x)*pi_prior_tmp+P_background_func(x)*(1-pi_prior_tmp))-postL_cutoff
    }
    root_tmp = uniroot(postL_func, lower = mu_0i[i], upper = 2*exp(alphai_beta_g[i]), postL_cutoff = .9)$root
    return(root_tmp)
  })


  return(matrix(count_cufoffs,nrow = nrow(alphai_beta_g)))
}


