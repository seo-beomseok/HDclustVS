#' Generating multi-dimensional Gaussian random partition.
#'
#' This function generates multi-dimensional Gaussian random partition with noisy variables for simulation study.
#' @param n Size of observations.
#' @param mean Mean vectors.
#' @param sigma Covariance matrix.
#' @param df Degree of freedom.
#' @param type Type of the noncentral multivariate t distribution. Refer to mvtnorm::rmvt().
#' @param lower Lower bound of truncated t-distribution.
#' @param upper Upper bound of truncated t-distribution.
#' @return truncated t-distribution samples
#' @importFrom mvtnorm rmvnorm
#' @export

########################################################
####  simulate data
########################################################

rtmvt <- function(n,mean,sigma,df=3,type='shifted',lower=rep(-inf,length(mean)),upper=rep(inf,length(mean))){
  
  truncate <- function(x,lower,upper){
    lowermat = matrix(lower,nrow(x),ncol(x))
    uppermat = matrix(upper,nrow(x),ncol(x))
    x = x[which(apply(x>lowermat & x<uppermat,1,prod)==1),]
    return(x)
  }
  
  x = rmvt(n, delta=mean, sigma=sigma, df=df, type=type)
  x = truncate(x,lower,upper)
  nsamp = nrow(x)
  iter = 1
  while(nsamp < n){
    xadd = rmvt(n-nsamp, delta=mean, sigma=sigma, df=df, type=type)
    xadd = truncate(xadd,lower,upper)
    x = rbind(x,xadd)
    nsamp = nrow(x)

    iter=iter+1
  }
  return(x)
}

