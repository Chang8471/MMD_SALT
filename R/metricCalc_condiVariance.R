

#' Calculate conditional variance
#'
#' @description Define variance metric as average squared distance from median to each point greater than the median
#'
#' @param X matrix of data, NA if not in the inclusion condition
#' @param margin which margin to apply function, 1 if rows are genes/probes
#'
#' @return marginal variance
#' @export
#'
#' @examples
applyVar = function(X,margin=1){
  tmp = apply(X,margin,function(x){
    x_median = median(x,na.rm = T)
    x_sq = (x-x_median)^2
    mean(x_sq[(x-x_median)>0],na.rm = T)
  })
  return(tmp)
}

#' Conditional variance pooled by group
#'
#' @description Define variance metric as average squared distance from median to each point greater than the median.
#'
#' @param X matrix of data, NA if not in the inclusion condition
#' @param g group labels, 2 levels
#' @param margin which margin to apply function, 1 if rows are genes/probes
#'
#' @return Conditional variance for each gene pooled by group
#' @export
#'
#' @examples
applyVar_groups = function(X,g,margin=1){
  tmp = apply(X,margin,function(x){
    ss = sum(tapply(x,g,function(y){
      y_median = median(y,na.rm = T)
      y_sq = (y-y_median)^2
      sum(y_sq[(y-y_median)>0],na.rm = T)
    }),na.rm = T)
    df = sum(tapply(x,g,function(y){
      y_median = median(y,na.rm = T)
      sum((y-y_median)>0,na.rm = T)
    }),na.rm = T)
    return(ss/df)
  })
  return(tmp)
}

#' Difference in group medians
#'
#' @param X matrix of data, NA if not in the inclusion condition
#' @param g group labels, 2 levels
#' @param margin which margin to apply function, 1 if rows are genes/probes
#'
#' @return Difference in group medians, 2nd level - 1st level
#' @export
#'
#' @examples
deltaMedian = function(X,g,margin = 1){
  tmp = apply(X,margin,function(x){
    m = tapply(x,g,median,na.rm=T)
    return(m[2]-m[1])
  })
  return(tmp)
}


#' T-like statistics
#'
#' @description calculated as change in group medians/SE, were SE is calculated with pooled conditional variance and df.
#'
#' @param X matrix of data, NA if not in the inclusion condition
#' @param g group labels, 2 levels
#' @param margin which margin to apply function, 1 if rows are genes/probes
#'
#' @return T-like statistics
#' @export
#'
#' @examples
tLikeStat = function(X, g, margin = 1){
  SE_tmp = tmp = apply(X,margin,function(x){
    ss = sum(tapply(x,g,function(y){
      y_median = median(y,na.rm = T)
      y_sq = (y-y_median)^2
      sum(y_sq[(y-y_median)>0],na.rm = T)
    }),na.rm = T)
    df = sum(tapply(x,g,function(y){
      y_median = median(y,na.rm = T)
      sum((y-y_median)>0,na.rm = T)
    }),na.rm = T)
    return(sqrt(ss/df)/sqrt(df))
  })
  deltaMedian = apply(X,margin,function(x){
    m = tapply(x,g,median,na.rm=T)
    return(m[2]-m[1])
  })
  return(deltaMedian/SE_tmp)

}
