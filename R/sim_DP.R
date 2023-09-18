#' Simulate nCounter samples and selectively introduce differential presence to subset of probes
#'
#' @param model_fit optional input, previously estimated MMD model parameters that contains study of interest
#' @param sampleGroup vector same length as total number of sampled used to fit ``model_fit``. For samples that's not in study of interest will have value NA, and the two groups in the study will have value 0 or 1 to indicate group assignment.
#' @param baseline_SigPrev baseline signal prevalence
#' @param effectSize target effect size for differential presence, will be added on top of ``baseline_SigPrev``
#' @param nPosControl number of positive control probes, default = 10
#' @param nNegControl number of negative control probes, default = 10
#' @param nSample total number of samples, default = 500
#' @param d0 prior degree of freedom for simulating signal variance
#' @param beta_offset a numerical value to offset all beta_offset, the larger the further away all signals will be from the noise
#' @param seed.par seed for simulating parameters, default = 0
#' @param seed.data seed for simulating Z and Y, default = 0
#' @param expressPercent fixed proportion of endogenous probes expressed in each sample
#'
#' @return
#' @export
#'
#'
#' @examples

sim_DP_wRealData = function(model_fit = NULL, sampleGroup = NULL,
                  baseline_SigPrev = .3, effectSize = .3,
                  nPosControl = 10, nNegControl = 10,
                  nSample = 500, d0 = 50, expressPercent = 1/5,beta_offset = 0,
                  seed.par=0, seed.data=0){

  set.seed(seed.par)


  #overwrite number of engougenous genes
  nEndogenous = length(model_fit$beta_g)
  # fix beta_i
  beta_g = c(model_fit$beta_g, 1:nNegControl, 1:nPosControl)
  # bootstrap alpha_i
  alpha_i = c(model_fit$alpha_i[sampleGroup==0&!is.na(sampleGroup)],model_fit$alpha_i[sampleGroup==1&!is.na(sampleGroup)],
              sample(model_fit$alpha_i, size = nSample-sum(!is.na(sampleGroup)),replace = T))


  nGene = nEndogenous+nPosControl+nNegControl


  # parameters for signal component
  sigma2_10 = exp(beta_g/1.3-4)
  sigma2_1g = (d0*sigma2_10)/rchisq(nGene,d0)
  sigma2_1g[(nGene-nPosControl+1):nGene]=0 # pos control
  alphai_beta_g = matrix(0,nrow = nGene,ncol = nSample)
  alphai_beta_g = sweep(alphai_beta_g,1,beta_g,"+")
  alphai_beta_g = sweep(alphai_beta_g,2,alpha_i,"+")

  # background parameters
  mu_0i = rlnorm(nSample, 3,.1)
  sigma_0i = mu_0i*exp(-.9)

  # Z_gi, fix a proportion of genes expressed in each sample
  Z_mtx = matrix(0,nGene,nSample)
  tmp = rep(0,nEndogenous)
  tmp[1:floor(expressPercent*nEndogenous)]=1
  Z_mtx[1:nEndogenous,1:nSample] = replicate(nSample,sample(tmp,length(tmp),replace = F)) # permute z label to decide which gene is detected in a sample
  Z_mtx[(nGene-nPosControl+1):nGene,] =1 # positive control
  Z_mtx[(nGene-nPosControl-nNegControl+1):(nGene-nPosControl),] =0 # neg control



  ########## simulate diff detection
  nGroupA = sum(sampleGroup==0,na.rm = T) # samples 1st group
  nGroupB = sum(sampleGroup==1,na.rm = T) # samples 2nd group
  nDDgenes = 100 # diffPres genes, half of total for symetric

  # get samples and genes of interest
  ind_groupA = 1:nGroupA
  ind_groupB = (1+nGroupA):(nGroupB+nGroupA)
  # sort genes by beta, select to add DD
  ind_genes_DD_1 = (rank(-beta_g[1:nEndogenous]) %%(nEndogenous/nDDgenes))==1
  ind_genes_DD_2 = (rank(-beta_g[1:nEndogenous]) %%(nEndogenous/nDDgenes))==2

  ########### up to here, all parameters simulated

  ########### simulate pi for each group
  set.seed(seed.data) # set seed for this run

  # simulate z matrix
  # group 1
  pi_tmp = matrix(baseline_SigPrev, nrow = nEndogenous, ncol=nGroupA)
  Z_mtx[(1:nEndogenous),ind_groupA] =
    matrix(rbinom(nEndogenous*nGroupA,1,
                  sweep(pi_tmp,1,ind_genes_DD_1*effectSize,"+")),
           nrow = nEndogenous, ncol=nGroupA)
  #group 2
  pi_tmp = matrix(baseline_SigPrev, nrow = nEndogenous, ncol=nGroupB)
  Z_mtx[(1:nEndogenous),ind_groupB] =
    matrix(rbinom(nEndogenous*nGroupB,1,
                  sweep(pi_tmp,1,ind_genes_DD_2*effectSize,"+")),
           nrow = nEndogenous, ncol=nGroupB)



  # simulate Y|Z
  Y_signal = matrix(rlnorm(nGene*nSample,c(alphai_beta_g),sqrt(sigma2_1g)),nrow = nGene)
  Y_background = matrix(rnorm(nSample*nGene,mu_0i,sigma_0i),nrow=nGene, byrow = T)
  Y_mtx_wCtrl = ifelse(Z_mtx, Y_signal+Y_background, Y_background)
  Y_mtx_wCtrl = ceiling(Y_mtx_wCtrl)
  Y_mtx_wCtrl[Y_mtx_wCtrl<1] = 1

  Y_mtx = Y_mtx_wCtrl[1:nEndogenous,]
  colnames(Y_mtx)=paste0("sample",1:nSample)
  rownames(Y_mtx)=paste0("gene",1:nEndogenous)

  data_tmp = list(
    Y_mtx = Y_mtx,
    Z_mtx = Z_mtx[1:nEndogenous,],
    Neg_control = Y_mtx_wCtrl[(nGene-nPosControl-nNegControl+1):(nGene-nPosControl),],
    Pos_control = Y_mtx_wCtrl[(nGene-nPosControl+1):nGene,],

    pi_i = expressPercent,mu_0i = mu_0i,sigma_0i=sigma_0i,
    alpha_i = alpha_i, beta_g = beta_g[1:nEndogenous],alpha_i_beta_g=alphai_beta_g[1:nEndogenous,],
    sigma2_10 = sigma2_10[1:nEndogenous], sigma2_1g = sigma2_1g[1:nEndogenous],

    ind_groupA=ind_groupA,ind_groupB=ind_groupB,
    ind_genes_DD_1=ind_genes_DD_1,ind_genes_DD_2=ind_genes_DD_2
  )


  return(data_tmp)


}
