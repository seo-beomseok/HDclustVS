#' Generating multi-dimensional Gaussian random partition.
#'
#' This function generates multi-dimensional Gaussian random partition with noisy variables for simulation study.
#' @param n Size of observations.
#' @param relev.blocks A vector of positive integers indicating the size of each relevant block.
#' @param irrel.blocks A vector of positive integers indicating the size of each irrelevant block.
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
#' sim3 = genData3(n=1000,relev.blocks=c(200,200,200,200,200),irrel.blocks=1000,relev.var=c(50,50),irrel.var=c(50,50),C=5,rep=1) 
#' @export

########################################################
####  simulate data
########################################################

genData3 <- function(n,relev.blocks,irrel.blocks,C,relev.var=c(100,200),irrel.var=c(100,200),df=4,overlap=FALSE,BClevel=10^-5,MinSize=1/(5*C),rep=1){
  while(1){
    pp = sort(runif(C-1,0,1))
    pp_true <- c(pp,1) - c(0,pp)   ## pi
    if(prod(pp_true>MinSize)) break
  }
  
  rsample <- rmultinom(n, 1, pp_true)      ## C by n, indicator of sample
  z_true <- apply(rsample, 2, function(x)which(x == 1))
  
  iter = 1
  p1 = sum(relev.blocks)
  
  p1h = round(p1/2,0)
  mu_relev <- matrix(0, nrow = p1, ncol = C)
  mu_relev[,1] = 4.0
  mu_relev[,2] = -4.0
  mu_relev[1:p1h,3] = -2; mu_relev[(p1h+1):p1,3] = 2
  mu_relev[1:p1h,4] = 2; mu_relev[(p1h+1):p1,4] = -2
  
  while(1){
    Sigma_relev = array(0, dim = c(p1,p1,C))    ## var-cov matrix of x,y coordinates 
    
    for(k in 1:length(relev.blocks)){
      r0 = ifelse(k==1,1,sum(relev.blocks[1:(k-1)])+1)
      r1 = sum(relev.blocks[1:k])
      for(i in 1:C){
        Sigma_relev[r0:r1,r0:r1,i] = genPositiveDefMat(dim=relev.blocks[k], covMethod = "onion", rangeVar = relev.var)$Sigma
      }
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
    if(!sum(bhattacharyyaCoef>BClevel,na.rm=T)) break
    if(iter>10 & !sum(bhattacharyyaCoef>10^-2,na.rm=T)){
      warning('In 10 iterations, generating non-overlapped clusters failed. The returned clusters are overlapped each other.')
      break
    }
    
    iter=iter+1
  }
  
  p2 = sum(irrel.blocks)
  X_irrel = matrix(0,n,p2)
  mu_irrel = runif(p2,-10,10)
  #Sigma_irrel = diag(runif(p2,irrel.var[1],irrel.var[2]))
  Sigma_irrel = genPositiveDefMat(dim=p2, covMethod = "onion", rangeVar = irrel.var)$Sigma
  
  ret = list()
  for(i in 1:rep){
    num_z = table(z_true)
    X_relev = matrix(NA,n,p1)
    for(k in 1:length(relev.blocks)){
      r0 = ifelse(k==1,1,sum(relev.blocks[1:(k-1)])+1)
      r1 = sum(relev.blocks[1:k])
      for(j in 1:C){
        X_relev[which(z_true==j),r0:r1] = rtmvt(n=num_z[j], mean=mu_relev[r0:r1,j], sigma=Sigma_relev[r0:r1,r0:r1,j], df=df, lower=rep(-100,relev.blocks[k]), upper=rep(100,relev.blocks[k]))
      }
    }
    
    X_irrel = rtmvt(n=n,mean=mu_irrel,sigma=Sigma_irrel, df=df, lower=rep(-100,p2), upper=rep(100,p2))
    X_total=cbind(X_relev,X_irrel)
    
    colnames(X_relev) = 1:p1
    colnames(X_irrel) = 1:p2
    colnames(X_total) = 1:(p1+p2)
    ret[[i]] = list(X_total=X_total,X_relev=X_relev,X_irrel=X_irrel,z=z_true)
  }  
  
  return(ret)
}  

