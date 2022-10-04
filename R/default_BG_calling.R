
#' Call background with default methods
#'
#' @description call background according to only negative control probes
#'
#' @param Ymtx Matrix of raw counts, gene by sample. Endogenous probes only
#' @param Neg_control Matrix of raw counts, gene by sample. Negative control probes only
#' @param method ways for calling background, "max" for less than max of all negative controls; "norm" model negative controls with normal distributions; "logNorm" model negative controls with log normal distributions
#' @param binary boolean indicator for returning binary calling (1 for signal and 0 for background). `method = "max"` will always return binary labels. If false, the other two method will return p-values for observing each count from sample specific background distribution.
#'
#' @return matrix of binary labels (1 for signal and 0 for background), or p-values for observing each count from sample specific background distribution. Same dimensions as ``Ymtx``
#' @export
#'
#' @examples
call_BG_wNegCtrl = function(Ymtx, Neg_control, method, binary){
  if (method == "max"){
    return(call_Z_maxNegCtrl(Ymtx, Neg_control))
  }else if (method == "norm"){
    Z_mtx = call_Z_normNegCtrl(Ymtx, Neg_control)
  }else if (method == "logNorm"){
    Z_mtx = call_Z_logNormNegCtrl(Ymtx, Neg_control)
  }
  if (binary) return(Z_mtx<.05)
  return(Z_mtx)
}

#' Call background with default methods, Z score
#'
#' @description call background according to only negative control probes
#'
#' @param Ymtx Matrix of raw counts, gene by sample. Endogenous probes only
#' @param Neg_control Matrix of raw counts, gene by sample. Negative control probes only
#' @param method ways for calling background, "max" for less than max of all negative controls; "norm" model negative controls with normal distributions; "logNorm" model negative controls with log normal distributions
#' @param binary cutoff used binarize Z score. `method = "max"` will always return binary labels. If false, the other two method will return p-values for observing each count from sample specific background distribution.
#'
#' @return matrix of binary labels (1 for signal and 0 for background), or Z scores for observing each count from sample specific background distribution. Same dimensions as ``Ymtx``
#' @export
#'
#' @examples
call_BG_wNegCtrl_Zscore = function(Ymtx, Neg_control, method, binary=NULL){
  if (method == "max"){
    return(call_Z_maxNegCtrl(Ymtx, Neg_control))
  }else if (method == "norm"){
    Z_mtx = call_Z_normNegCtrl(Ymtx, Neg_control)
  }else if (method == "logNorm"){
    Z_mtx = call_Z_logNormNegCtrl(Ymtx, Neg_control)
  }
  if (binary) return(Z_mtx>binary)
  else return(Z_mtx)
}



#' Call background with default methods
#'
#' @description call background according to only negative control probes
#'
#' @param Ymtx Matrix of raw counts, gene by sample. Endogenous probes only
#' @param Neg_control Matrix of raw counts, gene by sample. Negative control probes only
#'
#' @return matrix of binary labels (1 for signal and 0 for background) Same dimensions as ``Ymtx``
#' @noRd
call_Z_maxNegCtrl = function(Ymtx, Neg_control){
  Z_mtx = sapply(1:ncol(Ymtx),function(i){
    ifelse(Ymtx[,i]>max(Neg_control[,i]),1,0)
  })
  return(Z_mtx)
}

#' Call background with default methods
#'
#' @description call background according to only negative control probes
#'
#' @param Ymtx Matrix of raw counts, gene by sample. Endogenous probes only
#' @param Neg_control Matrix of raw counts, gene by sample. Negative control probes only
#' @param zScore boolean indicator to return Z score instead of p.value, default to F
#'
#' @return matrix of p-values for observing each count from sample specific background distribution. Same dimensions as ``Ymtx``
#' @noRd
call_Z_normNegCtrl = function(Ymtx, Neg_control, zScore = F){
  mu_tmp = colMeans(Neg_control)
  sd_tmp = apply(Neg_control,2,sd)
  if (zScore)  Z_mtx_est = sweep(sweep(Ymtx,2,mu_tmp,"-"),2,sd_tmp,"/")
  else  Z_mtx_est = t(pnorm(t(Ymtx),mu_tmp,sd_tmp,lower.tail = F))
  return(Z_mtx_est)
}

#' Call background with default methods
#'
#' @description call background according to only negative control probes
#'
#' @param Ymtx Matrix of raw counts, gene by sample. Endogenous probes only
#' @param Neg_control Matrix of raw counts, gene by sample. Negative control probes only
#' @param zScore boolean indicator to return Z score instead of p.value, default to F
#'
#' @return matrix of p-values for observing each count from sample specific background distribution. Same dimensions as ``Ymtx``
#' @noRd
call_Z_logNormNegCtrl = function(Ymtx, Neg_control, zScore = F){
  mu_tmp = colMeans(log(Neg_control))
  sd_tmp = apply(log(Neg_control),2,sd)

  if (zScore)  Z_mtx_est = sweep(sweep(log(Ymtx),2,mu_tmp,"-"),2,sd_tmp,"/")
  else  Z_mtx_est = t(pnorm(t(log(Ymtx)),mu_tmp,sd_tmp,lower.tail = F))

  return(Z_mtx_est)
}




