

#' Simulate NanoString nCounter data with from null model
#'
#' @description simulate parameter and data, with potentially varying prior probability of signal in all samples
#'
#' @param nEndogenous number of endogenous probes, default = 800
#' @param nPosControl number of positive control probes, default = 10
#' @param nNegControl number of negative control probes, default = 5
#' @param nSample total number of samples, default = 500
#' @param alwaysOnGenes for the endogenous probes, set certain amount always expressed in all samples, default = 0
#' @param fixedBackground boolean indicator, if true, all samples will share the same background mean and variance parameters
#' @param d0 prior degree of freedom for simulating signal variance
#' @param alpha_Offset a numerical value to offset all alpha_i, the larger the further away all signals will be from the noise
#'
#' @return a list of data and simulated true parameters
#'
#' \itemize{
#' \item `Y_mtx` a matrix of all observed counts for endogenous probes. Probes by samples
#' \item `Z_mtx` a matrix of true labels for endogenous probes, 0 for noise, 1 for signal, same dimensions as `Y_mtx`
#' \item `Neg_control`  a matrix of all observed counts for Negative control probes. Probes by sample
#' \item `Pos_control`  a matrix of all observed counts for Positive control probes. Probes by sample
#' \item `pi_i` vector of simulated sample marginal prior probabilty of observing signal
#' \item `mu_0i` vector simulated background mean
#' \item `sigma_0i` vector of simulate background variance
#' \item `alpha_i` vector of simulated sample contribution to signal mean, offset included
#' \item `beta_g` vector of simulated probe contribution to signal mean
#' \item `alpha_i_beta_g` matrix of signal mean, same dimensions as `Y_mtx`
#' \item `sigma2_10` vector of signal variance prior
#' \item `sigma2_1g` vector of simulated signal variance
#' }
#' @export
#'
#' @examples
simData = function(nEndogenous = 800, nPosControl = 10, nNegControl = 5,
                   nSample = 500, alwaysOnGenes = 0,
                   fixedBackground = T, d0 = 50, beta_offset = 0,
                   seed.par=NULL, seed.data=NULL){


  nGene = nEndogenous+nPosControl+nNegControl


  # parameters for signal component
  if(seed.par) set.seed(seed.par)
  alpha_i = rnorm(nSample,0,.1)
  beta_g = rnorm(nGene,3.4,2/3)+ beta_offset+mean(alpha_i)
  alpha_i = alpha_i-mean(alpha_i)
  beta_g[(nGene-nPosControl+1):nGene] = 1:nPosControl # pos control, gradient
  sigma2_10 = exp(beta_g/1.3-4)
  sigma2_1g = (d0*sigma2_10)/rchisq(nGene,d0)
  sigma2_1g[(nGene-nPosControl+1):nGene]=0 # pos control
  alphai_beta_g = matrix(0,nrow = nGene,ncol = nSample)
  alphai_beta_g = sweep(alphai_beta_g,1,beta_g,"+")
  alphai_beta_g = sweep(alphai_beta_g,2,alpha_i,"+")

  # background parameters
  mu_0i = rlnorm(nSample, 2.43,.44)
  if (fixedBackground) mu_0i = rep(20,nSample)
  sigma_0i = rnorm(nSample, mu_0i*exp(-.9),.03*mu_0i)

  # Z_gi, fix a proportion of genes expressed in each sample
  if(seed.data) set.seed(seed.data)
  Z_mtx = matrix(0,nGene,nSample)
  tmp = rep(0,nEndogenous-alwaysOnGenes)
  tmp[1:floor(expressPercent*(nEndogenous-alwaysOnGenes))]=1
  Z_mtx[1:(nEndogenous-alwaysOnGenes),1:nSample] = replicate(nSample,sample(tmp,length(tmp),replace = F))
  if (alwaysOnGenes>0){Z_mtx[(nEndogenous-alwaysOnGenes+1):nEndogenous,] =1} # always expressed in all samples, but have variance
  Z_mtx[(nGene-nPosControl+1):nGene,] =1 # positive control
  Z_mtx[(nGene-nPosControl-nNegControl+1):(nGene-nPosControl),] =0 # neg control


  # simulate Y|Z
  Y_signal = matrix(rlnorm(nGene*nSample,c(alphai_beta_g),sqrt(sigma2_1g)),nrow = nGene)
  Y_background = matrix(rnorm(nSample*nGene,mu_0i,sigma_0i),nrow=nGene, byrow = T)
  Y_mtx_wCtrl = ifelse(Z_mtx, Y_signal+Y_background, Y_background)
  Y_mtx_wCtrl = ceiling(Y_mtx_wCtrl)
  Y_mtx_wCtrl[Y_mtx_wCtrl<1] = 1

  Y_mtx = Y_mtx_wCtrl[1:nEndogenous,]
  colnames(Y_mtx)=paste0("sample",1:nSample)
  rownames(Y_mtx)=paste0("gene",1:nEndogenous)

  return(list(
    Y_mtx = Y_mtx,
    Z_mtx = Z_mtx[1:nEndogenous,],
    Neg_control = Y_mtx_wCtrl[(nGene-nPosControl-nNegControl+1):(nGene-nPosControl),],
    Pos_control = Y_mtx_wCtrl[(nGene-nPosControl+1):nGene,],

    pi_i = pi_i,mu_0i = mu_0i,sigma_0i=sigma_0i,
    alpha_i = alpha_i, beta_g = beta_g[1:nEndogenous],alpha_i_beta_g=alphai_beta_g[1:nEndogenous,],
    sigma2_10 = sigma2_10[1:nEndogenous], sigma2_1g = sigma2_1g[1:nEndogenous]
  ))
}



#' Simulate from null model, with fixed proportion of endogenous probes having signal
#'
#' @description Simulate parameters and data. While all samples have same proportion of probes having signal, it is a different subset for each sample.
#'
#' @param nEndogenous number of endogenous probes, default = 800
#' @param nPosControl number of positive control probes, default = 10
#' @param nNegControl number of negative control probes, default = 5
#' @param nSample total number of samples, default = 500
#' @param alwaysOnGenes for the endogenous probes, set certain amount always expressed in all samples, default = 0
#' @param fixedBackground boolean indicator, if true, all samples will share the same background mean and variance parameters
#' @param d0 prior degree of freedom for simulating signal variance
#' @param alpha_Offset a numerical value to offset all alpha_i, the larger the further away all signals will be from the noise
#' @param seed.ab seed for simulating parameters for signal components
#' @param seed.sim seed for simulating background parameters, Z and Y
#' @param expressPercent fixed proportion of endogenous probes expressed in each sample
#' @param alpha_shape1 first parameter to simulated alpha_i from rbeta()
#' @param alpha_shape2 second parameter to simulated alpha_i from rbeta()
#'
#' @return a list of data and simulated true parameters
#'
#' \itemize{
#' \item `Y_mtx` a matrix of all observed counts for endogenous probes. Probes by samples
#' \item `Z_mtx` a matrix of true labels for endogenous probes, 0 for noise, 1 for signal, same dimensions as `Y_mtx`
#' \item `Neg_control`  a matrix of all observed counts for Negative control probes. Probes by sample
#' \item `Pos_control`  a matrix of all observed counts for Positive control probes. Probes by sample
#' \item `pi_i` true prior probability of signal in a sample, same as `expressPercent` input parameter
#' \item `mu_0i` vector simulated background mean
#' \item `sigma_0i` vector of simulate background variance
#' \item `alpha_i` vector of simulated sample contribution to signal mean, offset included
#' \item `beta_g` vector of simulated probe contribution to signal mean
#' \item `alpha_i_beta_g` matrix of signal mean, same dimensions as `Y_mtx`
#' \item `sigma2_10` vector of signal variance prior
#' \item `sigma2_1g` vector of simulated signal variance
#' }
#' @export
#'
#' @examples
simData_fixSigProp = function(nEndogenous = 800,
                              nPosControl = 10, nNegControl = 10,
                              nSample = 500, alwaysOnGenes = 0,
                              fixedBackground = F, d0 = 50, alpha_Offset = 0,
                              seed.ab = 0,seed.sim = 1, expressPercent = 1/5,
                              alpha_shape1=3,alpha_shape2 =3){



  nGene = nEndogenous+nPosControl+nNegControl


  # parameters for signal component
  set.seed(seed.ab)
  alpha_i = rbeta(nSample,alpha_shape1,alpha_shape2) +alpha_Offset
  beta_g = rnorm(nGene,3.4,2/3)
  beta_g[(nGene-nPosControl+1):nGene] = 1:nPosControl # pos control, gradient
  sigma2_10 = exp(beta_g/1.3-4)
  sigma2_1g = (d0*sigma2_10)/rchisq(nGene,d0)
  sigma2_1g[(nGene-nPosControl+1):nGene]=0 # pos control
  alphai_beta_g = matrix(0,nrow = nGene,ncol = nSample)
  alphai_beta_g = sweep(alphai_beta_g,1,beta_g,"+")
  alphai_beta_g = sweep(alphai_beta_g,2,alpha_i,"+")


  set.seed(seed.sim)
  # background parameters
  mu_0i = rlnorm(nSample, 3,.1)
  if (fixedBackground) mu_0i = rep(20,nSample)
  sigma_0i = mu_0i*exp(-.9)


  # Z_gi, fix a proportion of genes expressed in each sample
  Z_mtx = matrix(0,nGene,nSample)
  tmp = rep(0,nEndogenous-alwaysOnGenes)
  tmp[1:floor(expressPercent*(nEndogenous-alwaysOnGenes))]=1
  Z_mtx[1:(nEndogenous-alwaysOnGenes),1:nSample] = replicate(nSample,sample(tmp,length(tmp),replace = F))
  if (alwaysOnGenes>0){Z_mtx[(nEndogenous-alwaysOnGenes+1):nEndogenous,] =1} # always expressed in all samples, but have variance
  Z_mtx[(nGene-nPosControl+1):nGene,] =1 # positive control
  Z_mtx[(nGene-nPosControl-nNegControl+1):(nGene-nPosControl),] =0 # neg control


  # simulate Y|Z
  Y_signal = matrix(rlnorm(nGene*nSample,c(alphai_beta_g),sqrt(sigma2_1g)),nrow = nGene)
  Y_background = matrix(rnorm(nSample*nGene,mu_0i,sigma_0i),nrow=nGene, byrow = T)
  Y_mtx_wCtrl = ifelse(Z_mtx, Y_signal+Y_background, Y_background)
  Y_mtx_wCtrl = ceiling(Y_mtx_wCtrl)
  Y_mtx_wCtrl[Y_mtx_wCtrl<1] = 1

  Y_mtx = Y_mtx_wCtrl[1:nEndogenous,]
  colnames(Y_mtx)=paste0("sample",1:nSample)
  rownames(Y_mtx)=paste0("gene",1:nEndogenous)

  return(list(
    Y_mtx = Y_mtx,
    Z_mtx = Z_mtx[1:nEndogenous,],
    Neg_control = Y_mtx_wCtrl[(nGene-nPosControl-nNegControl+1):(nGene-nPosControl),],
    Pos_control = Y_mtx_wCtrl[(nGene-nPosControl+1):nGene,],

    pi_i = expressPercent,mu_0i = mu_0i,sigma_0i=sigma_0i,
    alpha_i = alpha_i, beta_g = beta_g[1:nEndogenous],alpha_i_beta_g=alphai_beta_g[1:nEndogenous,],
    sigma2_10 = sigma2_10[1:nEndogenous], sigma2_1g = sigma2_1g[1:nEndogenous]
  ))
}

