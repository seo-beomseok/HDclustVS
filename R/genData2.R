#' Generating multi-dimensional Gaussian partition based on given parameters.
#'
#' This function generates multi-dimensional Gaussian partition with noisy variables for simulation study based on pre-specified mean and covariance matrix parameters.
#' @param n Size of observations.
#' @param p1 Dimension of relevant variables.
#' @param p2 Dimension of irrelevant variables.
#' @param C The number of clusters.
#' @param rep The number of data sets. 
#' @return A list of the synthesized data matrix.
#' @importFrom clusterGeneration genPositiveDefMat
#' @importFrom mvtnorm rmvnorm
#' @examples 
#' sim1 = genData(n=300,p1=100,p2=100,C=5,rep=1) 
#' sim2 = genData2(n=300,p1=100,p2=100,C=5,rep=1) 
#' @export

########################################################
####  simulate data
########################################################

genData2 <- function(n,p1,p2,C,rep=10){
  pp_true = c(0.2,0.5,0.05,0.15,0.1)
  
  rsample <- rmultinom(n, 1, pp_true)
  z_true <- apply(rsample, 2, function(x)which(x == 1))
  
  mu_relev = matrix(NA,p1,C)
  mu_relev[,1] = rep(3,100)
  mu_relev[,2] = rep(-3,100)
  mu_relev[,3] = c(rep(-3,50),rep(3,50))
  mu_relev[,4] = c(rep(3,50),rep(-3,50))
  mu_relev[,5] = rep(0,100)
  
  Sigma_relev = array(0, dim = c(p1,p1,C))    ## var-cov matrix of x,y coordinates 
  
  {
    library(matrixcalc)
    Sigma_relev[,,1] = 0.9
    diag(Sigma_relev[,,1]) = 2
    
    print(is.positive.definite(Sigma_relev[,,1]))
    
    Sigma_relev[,,2] = 0.6
    Sigma_relev[1:2,,2] = -0.6
    Sigma_relev[,1:2,2] = -0.6
    diag(Sigma_relev[,,2]) = 2
    
    print(is.positive.definite(Sigma_relev[,,2]))
    
    Sigma_relev[,,3] = 0.2
    diag(Sigma_relev[,,3]) = 2
    
    print(is.positive.definite(Sigma_relev[,,3]))
    
    Sigma_relev[,,4] = 0.6
    Sigma_relev[1:2,,4] = -0.3
    Sigma_relev[,1:2,4] = -0.3
    diag(Sigma_relev[,,4]) = 3
    
    print(is.positive.definite(Sigma_relev[,,4]))
    
    Sigma_relev[,,5] = 0.6
    Sigma_relev[1:2,,5] = -0.3
    Sigma_relev[,1:2,5] = -0.3
    diag(Sigma_relev[,,5]) = 3
    
    print(is.positive.definite(Sigma_relev[,,5]))
  }
  
  X_irrel = matrix(0,n,p2)
  mu_irrel = runif(p2,-3,3)
  Sigma_irrel = diag(runif(p2,10,50))
  
  ret = list()
  for(i in 1:rep){
    X_relev <- t(sapply(z_true, function(x) {rmvnorm(1, mu_relev[,x], Sigma_relev[,,x])}))
    X_irrel = rmvnorm(n,mu_irrel,Sigma_irrel)
    X_total=cbind(X_relev,X_irrel)
    
    colnames(X_relev) = 1:p1
    colnames(X_irrel) = 1:p2
    colnames(X_total) = 1:(p1+p2)
    ret[[i]] = list(X_total=X_total,X_relev=X_relev,X_irrel=X_irrel,z=z_true)
  }  
  
  return(ret)
}  
