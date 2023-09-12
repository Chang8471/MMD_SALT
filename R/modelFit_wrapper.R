#' Wrapper function to fit GMM
#'
#' @description This function calls ``parEst_itr()`` and ``forceMono()``. It fits the MMD iteratively and by the end calculate a posterior likelihood for each observed count containing signal.
#'
#' @param Y_mtx Matrix of raw counts, gene by sample
#' @param NegCtrl Matrix of raw counts for negative control probes, gene by sample. Must have the same number of columns as ``Y_mtx``. Will be ignored if ``Z_mtx`` is provided
#' @param Z_mtx Optional matrix of binary labels, same dimensions as ``Y_mtx``. 0 as in background, and 1 for signal.
#' @param maxIter Maximum iterations to fit the model
#' @param verbose Print steps of model fitting if True
#' @param ... The other provided arguments will go into ``parEst_itr()``
#'
#'
#'
#' @return a list of objects
#' \itemize{
#' \item `model_fit` fitted MMD model. Output from ``parEst_itr()``
#' \item `Z_mtx_adj` Matrix of posterior likelihood after monotonicity adjustment. Ouput from ``forceMono()``
#' }
#' @export
#'
#' @examples

fit_MMD = function(Y_mtx, maxIter = 20, NegCtrl = NULL, Z_mtx = NULL, verbose = F, ...){

  # check required input
  if (is.null(NegCtrl) & is.null(Z_mtx)) stop("One of NegCtrl or Z_mtx need to be provided")
  else if (is.null(Z_mtx)) Z_mtx = call_BG_wNegCtrl_Zscore(Y_mtx,NegCtrl,method="logNorm",binary = 2) # initialize Z_mtx according to the negative controls


  # parameter estimation
  if (verbose) print("Fitting iteration 1")
  model_fit = parEst_itr(Y_mtx,Z_mtx=Z_mtx, ...)

  for(i in 2:20){
    if (verbose) print(paste("iter",i))
    model_fit_tmp = parEst_itr(Y_mtx,pi_mtx=model_fit[[1]],
                         EBrobust=F,filterBackground  =T,
                         Beta_bayes = F, Beta_0_weight = F, Beta_kappa = 10)
    if(identical((model_fit_tmp[[1]]<.5 ),(model_fit[[1]]<.5))){
      if (verbose) print(paste("Stopped after iteration: ",i))
      model_fit = model_fit_tmp
      break
    }
    model_fit = model_fit_tmp
  }


  if (verbose) print("Forcing monotonicity...")
  Z_mtx_adj = forceMono(Y_mtx,
                        useparEst_itr_out = F,
                        alphai_beta_g = model_fit$m_hat,
                        sigma2_1g = model_fit$S2_hat,
                        mu_0i = model_fit$mu_0i,
                        sigma2_0i = model_fit$sigma2_0i,
                        pi_prior = model_fit$pi_prior)

  if (verbose) print("Model fitting complete")

  return(list(model_fit = model_fit, Z_mtx_adj = Z_mtx_adj))
}
