#' calculate test statistics for LRT
#'
#' @param Z0
#' @param Z1
#' @param n0
#' @param n1
#' @param s0
#' @param s1
#'
#' @return
#' @export
#'
#' @examples
calc_S_LR  = function(Z0,Z1,n0,n1,s0,s1){
  # maximize log likelihood under H null
  tmp = optimize(function(p) sum(Z0*log(p*s0))+sum((1-Z0)*log(1-p*s0))+
                   sum(Z1*log(p*s1))+sum((1-Z1)*log(1-p*s1)), c(0, 1),
                 maximum = TRUE)
  p_hat = tmp$maximum; l0 = tmp$objective

  # maximize log likelihood under H alter
  tmp = optimize(function(p) sum(Z0*log(p*s0))+sum((1-Z0)*log(1-p*s0)), c(0, 1), maximum = TRUE)
  p0_hat = tmp$maximum; l1 = tmp$objective

  tmp = optimize(function(p) sum(Z1*log(p*s1))+sum((1-Z1)*log(1-p*s1)), c(0, 1), maximum = TRUE)
  p1_hat = tmp$maximum; l1 = l1 +tmp$objective

  # LR statistics for one-sided test
  return(c(LR_stat = -2*(l0-l1), p_diff = p1_hat-p0_hat,p0_hat = p0_hat, p1_hat = p1_hat))
}



#' p null MLE, using the known sensitivities
#'
#' @param Z0
#' @param Z1
#' @param n0
#' @param n1
#' @param s0
#' @param s1
#'
#' @return
#' @export
#'
#' @examples
p_null_MLE = function(Z0,Z1,n0,n1,s0,s1){
  # maximize log likelihood under H null
  tmp = optimize(function(p) sum(Z0*log(p*s0))+sum((1-Z0)*log(1-p*s0))+
                   sum(Z1*log(p*s1))+sum((1-Z1)*log(1-p*s1)), c(0, 1),
                 maximum = TRUE)
  return(tmp$maximum)
}


#' simulate null distribution for LRT statistics
#'
#' @param n0
#' @param n1
#' @param p0
#' @param s0
#' @param s1
#' @param nRep
#'
#' @return
#' @export
#'
#' @examples
S_LRT_null = function( n0,n1,p0,s0,s1,nRep = 1000){
  null_S = sapply(1:nRep, function(seed){
    set.seed(seed)
    # simulate data under null
    Z0 = rbinom(n0,1,p0*s0); Z1 = rbinom(n1,1,p0*s1)
    return(calc_S_LR(Z0,Z1,n0,n1,s0,s1))
  })
}
