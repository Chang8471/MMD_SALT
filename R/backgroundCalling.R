
parEst_itr = function(Y_mtx,pi_mtx = NULL,Z_mtx = NULL,
                      EBrobust=F,returnPar = F,adj_posterior=F,filterBackground  =F,
                      fixPiProportion = F,approxSignalPar = F,subtractBackground = F,
                      Beta_bayes = F, Beta_0_weight = F, Beta_kappa = 10,
                      shrinktoMean=F){
  # use the initial marginal pi_i, and do not update in iterations
  if(fixPiProportion){
    pi_i = colMeans(Z_mtx)
    # but rough call is still made according to last round's posterior Z
    if (!is.null(pi_mtx)) Z_mtx = pi_mtx>.5
  }else{# update pi_i, each iteration
    # initialize and estimate pi_i
    if (!is.null(pi_mtx)){
      Z_mtx = pi_mtx>.5
      pi_i = colMeans(pi_mtx)
    }else if (!is.null(Z_mtx)){# only initial calls is provided
      pi_i = colMeans(Z_mtx)
    }
  }

  # check if all elements are provided
  if (is.null(Z_mtx)&is.null(pi_mtx)){ message("error: Z or pi mtx not correctly provided");return(NULL)}

  ############## background
  Y_mtx_background = Y_mtx
  Y_mtx_background[Z_mtx==1] = NA
  mu_0i = matrix(colMeans(Y_mtx_background,T),ncol = ncol(Y_mtx),nrow = nrow(Y_mtx),byrow = T)
  sigma2_0i = matrix(apply(Y_mtx_background,2,var,na.rm=T),ncol = ncol(Y_mtx),nrow = nrow(Y_mtx),byrow = T)


  ################### take above background
  if(subtractBackground) {
    logY_mtx = log(Y_mtx-mu_0i)
    if(any(is.na(logY_mtx)&Z_mtx)) message("error: background subtraction produce NA after log transformation")
  }else  {logY_mtx = log(Y_mtx)}
  logY_mtx[Z_mtx==0] = NA

  x_gene = factor(rep(rownames(logY_mtx),ncol(logY_mtx)),
                  levels=rownames(logY_mtx))
  x_sample = factor(rep(colnames(logY_mtx),each=nrow(logY_mtx)),
                    levels=colnames(logY_mtx))
  logY_vec = c(logY_mtx)

  # make indicator x and y
  x_gene = x_gene[!is.na(logY_vec)]
  x_sample = x_sample[!is.na(logY_vec)]
  logY_vec = logY_vec[!is.na(logY_vec)]
  x_mtx = model.matrix(~x_gene+x_sample-1) # model matrix
  #dim(x_mtx) # 83666  1469

  # fit lm
  lm_fitted_tmp = biglm(logY_vec~x_gene+x_sample-1,
                        data = data.frame(logY_vec,x_gene,x_sample))

  beta = coef(lm_fitted_tmp) # estimated beta
  if(shrinktoMean) {beta[is.na(beta)] = mean(beta[1:nrow(Y_mtx)],na.rm = T)
  }else   beta[is.na(beta)] = min(beta[1:nrow(Y_mtx)],na.rm = T)####genes never above background???

  logY_fitted = x_mtx %*% cbind(beta) # fitted

  # above background mean matrix
  alphai_beta_g = matrix(0,nrow = nrow(Y_mtx),ncol=ncol(Y_mtx))
  alphai_beta_g = sweep(alphai_beta_g,1,beta[1:nrow(Y_mtx)],"+")
  alphai_beta_g = sweep(alphai_beta_g,2,c(0,beta[(1+nrow(Y_mtx)):length(beta)]),"+")

  # residual for each observation
  resid = logY_vec - logY_fitted

  # gene specific variance, EB estimation
  s2_g = tapply(resid^2, x_gene, mean)
  df_g = tapply(resid, x_gene, length)
  df_g[is.na(df_g)]=0 # not detected genes

  var_eBayes = squeezeVar(s2_g,df_g, beta[1:nrow(Y_mtx)], robust = EBrobust)
  sigma2_1g = matrix(var_eBayes$var.post,ncol = ncol(Y_mtx),nrow = nrow(Y_mtx))

  # add bayes prior on beta_g
  if(Beta_bayes){
    # estimate beta_g
    beta_tmp = coef(lm_fitted_tmp)[1:nrow(Y_mtx)]
    beta_0 = mean(beta_tmp,na.rm = T)
    if(Beta_0_weight) beta_0= sum(beta_tmp*df_g,na.rm = T)/sum(df_g)
    beta_g_post = (Beta_kappa*beta_0+df_g*beta[1:nrow(Y_mtx)])/(Beta_kappa+df_g)

    # estimate sigma^2_g
    tmpsq = (beta_tmp-beta_0)^2
    tmpsq[is.na(tmpsq)] = 0
    sigma2_1g1 = sigma2_1g[,1]+(Beta_kappa*df_g/(df_g+Beta_kappa)*tmpsq)/(var_eBayes$df.prior+df_g)


    beta[1:nrow(Y_mtx)] = beta_g_post
    # update Expectation matrix
    alphai_beta_g = matrix(0,nrow = nrow(Y_mtx),ncol=ncol(Y_mtx))
    alphai_beta_g = sweep(alphai_beta_g,1,beta[1:nrow(Y_mtx)],"+")
    alphai_beta_g = sweep(alphai_beta_g,2,c(0,beta[(1+nrow(Y_mtx)):length(beta)]),"+")
    # update variance matrix
    sigma2_1g = matrix(sigma2_1g1,ncol = ncol(Y_mtx),nrow = nrow(Y_mtx))

  }
  ############ posterior f(z=1|Y)
  P_signal = dlnorm(Y_mtx,alphai_beta_g,sqrt(sigma2_1g))
  m_hat=NULL
  S2_hat = NULL
  if (approxSignalPar){
    EY = mu_0i +exp(alphai_beta_g+sigma2_1g/2)
    VarY = sigma2_0i+(exp(sigma2_1g)-1)*exp(2*alphai_beta_g+sigma2_1g)
    S2_hat = log(VarY/(EY)^2+1)
    m_hat = log(EY)-S2_hat/2
    P_signal = dlnorm(Y_mtx,m_hat,sqrt(S2_hat))
  }
  P_background = dnorm(Y_mtx,mu_0i,sqrt(sigma2_0i))

  pi_prior = matrix(pi_i,ncol = ncol(Y_mtx), nrow = nrow(Y_mtx), byrow = T)
  pi_mtx_posterior = P_signal*pi_prior/(P_signal*pi_prior+P_background*(1-pi_prior))


  if (adj_posterior) {
    # set upper bound for Z_hat according to alpha_i+beta_g
    P_signal_bound = dlnorm(exp(alphai_beta_g),alphai_beta_g,sqrt(sigma2_1g))
    P_background_bound = dnorm(exp(alphai_beta_g),mu_0i,sqrt(sigma2_0i))
    pi_mtx_posterior_bound = P_signal_bound*pi_prior/(P_signal_bound*pi_prior+P_background_bound*(1-pi_prior))

    pi_mtx_posterior = ifelse(log(Y_mtx)<alphai_beta_g,pi_mtx_posterior,pi_mtx_posterior_bound)
  }
  if(filterBackground){
    # any count less than background mean +1.96*sigma is background
    pi_mtx_posterior = ifelse(Y_mtx<mu_0i+1.96*sqrt(sigma2_0i),0,pi_mtx_posterior)
  }

  # return
  if(returnPar) return(list(
    pi_mtx_posterior=pi_mtx_posterior, beta=beta,
    alphai_beta_g=alphai_beta_g,sigma2_1g=sigma2_1g,mu_0i=mu_0i,sigma2_0i=sigma2_0i,
    var_eBayes=var_eBayes,s2_g=s2_g,df_g=df_g,
    pi_prior = pi_prior,m_hat=m_hat,S2_hat = S2_hat))
  return(list(pi_mtx_posterior=pi_mtx_posterior, beta=beta))

}

forceMono = function(Y_mtx,useparEst_itr_out,  parEst_itr_out=NULL,gridResolution=500,
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
  # seq of nodes
  Y_nodes = exp(seq(0,log(max(Y_mtx)+1),length.out = gridResolution))

  P_signal_tmp = dlnorm(Y_nodes,
                        rep(alphai_beta_g,each = length(Y_nodes)),
                        rep(sqrt(sigma2_1g),each = length(Y_nodes)))
  P_background_tmp = dnorm(Y_nodes,
                           rep(mu_0i,each = length(Y_nodes)),
                           rep(sqrt(sigma2_0i),each = length(Y_nodes)))
  pi_prior_tmp = rep(pi_prior, each =length(Y_nodes))
  pi_mtx_posterior_tmp = P_signal_tmp*pi_prior_tmp/(P_signal_tmp*pi_prior_tmp+P_background_tmp*(1-pi_prior_tmp))
  #sum(is.na(pi_mtx_posterior_tmp)) # turn NA to 0, since we will cummax later
  pi_mtx_posterior_tmp[is.na(pi_mtx_posterior_tmp)] = 0
  pi_mtx_posterior_tmp = matrix(pi_mtx_posterior_tmp,nrow = length(Y_nodes))
  #dim(pi_mtx_posterior_tmp)

  Z_post_mono = apply(pi_mtx_posterior_tmp,2,cummax) # find z that are monotonic
  #dim(Z_post_mono)

  end1 = findInterval(c(Y_mtx), Y_nodes) # indices of two end of the interval in
  end2 = end1+1; end2[end2>length(Y_nodes)]=length(Y_nodes)
  z1= Z_post_mono[matrix(c(end1,1:length(end1)),ncol = 2)]
  z2 = Z_post_mono[matrix(c(end2,1:length(end1)),ncol = 2)]
  x1= log(Y_nodes[end1])
  x2= log(Y_nodes[end2])
  z_mono = z2+(z1-z2)/(x1-x2)*(c(log(Y_mtx)-x2))
  z_final = matrix(z_mono,nrow = nrow(Y_mtx))

  return(z_final)
}





#' Inverse logit
#'
#' @param x real number value
#'
#' @return 1/(exp(-x)+1)
#' @export
#'
#' @examples
inv.logit = function(x){1/(exp(-x)+1)}


#' Logit
#'
#' @param x value between 0 and 1
#'
#' @return log(x/(1-x))
#' @export
#'
#' @examples
logit = function(x){log(x/(1-x))}
