

simData = function(nEndogenous = 800, nPosControl = 10, nNegControl = 5,
                   nSample = 500, alwaysOnGenes = 0,
                   fixedBackground = T, d0 = 50, alpha_Offset = 2){

  inv.logit = function(x){1/(exp(-x)+1)}
  nGene = nEndogenous+nPosControl+nNegControl
  # background
  mu_0i = rlnorm(nSample, 3,.1)
  if (fixedBackground) mu_0i = rep(20,nSample)
  sigma_0i = mu_0i*exp(-.9)
  Y_background = matrix(rnorm(nSample*nGene,mu_0i,sigma_0i),nrow=nGene, byrow = T)

  # signal
  alpha_i = alpha_Offset+rlnorm(nSample,1,.1) # assume independence
  beta_g = rnorm(nGene,3.4,2/3)
  beta_g[(nGene-nPosControl+1):nGene] = 1:nPosControl # pos control, gradient
  sigma2_10 = exp(beta_g/1.3-4)
  sigma2_1g = (d0*sigma2_10)/rchisq(nGene,d0)
  sigma2_1g[(nGene-nPosControl+1):nGene]=0 # pos control
  alphai_beta_g = matrix(0,nrow = nGene,ncol = nSample)
  alphai_beta_g = sweep(alphai_beta_g,1,beta_g,"+")
  alphai_beta_g = sweep(alphai_beta_g,2,alpha_i,"+")
  Y_signal = matrix(rlnorm(nGene*nSample,c(alphai_beta_g),sqrt(sigma2_1g)),nrow = nGene)

  # Z_gi
  pi_i = inv.logit(rnorm(nSample,-1.3,.4))
  pi_i[pi_i<0] = 0; pi_i[pi_i>1] = 1
  Z_mtx = matrix(rbinom(nSample*nGene,1,pi_i), nrow = nGene,byrow = T)
  if (alwaysOnGenes>0){Z_mtx[(nEndogenous-alwaysOnGenes+1):nEndogenous,] =1} # always expressed in all samples, but have variance
  Z_mtx[(nGene-nPosControl+1):nGene,] =1 # positive control
  Z_mtx[(nGene-nPosControl-nNegControl+1):(nGene-nPosControl),] =0 # neg control

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


simData_fixSigProp = function(nEndogenous = 800,
                              nPosControl = 10, nNegControl = 10,
                              nSample = 500, alwaysOnGenes = 0,
                              fixedBackground = F, d0 = 50, alpha_Offset = 0,
                              seed.ab = 0,seed.sim = 1, expressPercent = 1/5,
                              alpha_shape1=3,alpha_shape2 =3){




  inv.logit = function(x){1/(exp(-x)+1)}
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

