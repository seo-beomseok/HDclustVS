#' Generating multi-dimensional Gaussian random partition.
#'
#' This function generates multi-dimensional Gaussian random partition with noisy variables for simulation study.
#' @param n Size of observations.
#' @param p1 Dimension of relevant variables.
#' @param p2 Dimension of irrelevant variables.
#' @param C The number of clusters.
#' @param relev.var Range of positive number. Default is [10,50]. Range for variances of relevant variables.
#' @param irrel.var Range of positive number. Default is [100,400]. Range for variances of irrelevant variables.
#' @param overlap Logical. If it is True, the algorithm allows clusters on relevant variable space to overlap each other. Overlap is evaluated by Bhattacharyya coefficient.
#' @param BClevel Threshold of Bhattacharyya coefficient. Default is 10^-5.
#' @param MinSize Minimum proportion allowed for each cluster. Default is 1/(5*C).
#' @return A synthesized data matrix.
#' @importFrom clusterGeneration genPositiveDefMat
#' @importFrom mvtnorm rmvnorm
#' @examples 
#' sim1 = genData(n=300,p1=100,p2=100,C=5,rep=1)
#' sim2 = genData2(n=300,p1=100,p2=100,C=5,rep=1) 
#' @export

########################################################
####  simulate data
########################################################

genData <- function(n,p1,p2,C,relev.var=c(10,50),irrel.var=c(100,400),overlap=FALSE,BClevel=10^-5,MinSize=1/(5*C),rep=1){
  while(1){
    pp = sort(runif(C-1,0,1))
    pp_true <- c(pp,1) - c(0,pp)   ## pi
    if(prod(pp_true>MinSize)) break
  }
  
  rsample <- rmultinom(n, 1, pp_true)      ## C by n, indicator of sample
  z_true <- apply(rsample, 2, function(x)which(x == 1))
  
  while(1){
    mu_relev <- matrix(runif(p1*C,-10,10), nrow = p1, ncol = C)
    Sigma_relev = array(0, dim = c(p1,p1,C))    ## var-cov matrix of x,y coordinates 
    for(i in 1:C){
      Sigma_relev[,,i] = genPositiveDefMat(dim=p1, covMethod = "onion", rangeVar = relev.var)$Sigma
    }
    
    bhattacharyyaDist = matrix(NA,C,C)
    for(i in 1:(C-1)){
      for(j in (i+1):C){
        mu1 = mu_relev[,i]
        mu2 = mu_relev[,j]
        dif = mu1-mu2
        sig1 = Sigma_relev[,,i]
        sig2 = Sigma_relev[,,j]
        sig.avg = (sig1+sig2)/2
        sig.inv = solve(sig.avg)
        bhattacharyyaDist[i,j] = (1/8)*(t(dif)%*%sig.inv%*%dif)+(1/2)*log(det(sig.avg)/sqrt(det(sig1)*det(sig2)))
      }
    }
    bhattacharyyaCoef = exp(-bhattacharyyaDist)
    if(overlap==TRUE) break
    if(!sum(bhattacharyyaCoef>10^-5,na.rm=T)) break
  }
  
  X_irrel = matrix(0,n,p2)
  mu_irrel = runif(p2,-10,10)
  Sigma_irrel = diag(runif(p2,irrel.var[1],irrel.var[2]))
  
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

