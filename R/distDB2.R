#' # Davis-Bouldin Type Distance
#'
#' This function calculates Davis-Bouldin type distance from the estimated GMM components of the output of HMM-VB.
#' @param X A data matrix.
#' @param hmmvb An object of class 'HMMVB'
#' @param clust An object of class 'HMMVBclust'
#' @return A distance matrix.
#' @export

########################################################
####  Davis-Bouldin type distance
########################################################

distDB2 = function(X,hmmvb,clust){
  KK = hmmvb@VbStructure@numst
  TT = hmmvb@VbStructure@nb
  NN = dim(X)[1]
  BD = hmmvb@VbStructure@bdim     # size of each block
  HC = hmmvb@HmmChain             # hmm chain
  VV = t(matrix(unlist(clust@clustParam$vseq[clust@clsid]),TT))+1
  
  for(t in 1:TT){
    if(KK[t]==1){
      HC[[t]]@mean = t(matrix(rep(HC[[t]]@mean,max(KK)),BD[t],max(KK)))
      HC[[t]]@sigma = lapply(1:max(KK),function(x) HC[[t]]@sigma[[1]])
    }
  }
  
  WS = array(0,dim=c(max(KK),max(KK),TT))
  for(t in 1:TT){
    MW = HC[[t]]@mean
    if(BD[t]==1){
      MS = t(t(sapply(HC[[t]]@sigma,diag)))
    }else{
      MS = t(sapply(HC[[t]]@sigma,diag))
    }
    
    for(i in 1:(KK[t]-1)){
      for(j in (i+1):KK[t]){
        WS[j,i,t] = WS[i,j,t] = sqrt(sum(((MW[i,]-MW[j,])^2)/(MS[i,]+MS[j,])))
      }
    }
    
  }
  
  DD = rep(NA,NN*(NN-1)/2)
  for(i in 1:(NN-1)){
    for(j in (i+1):NN){
      DD[NN*(i-1)-i*(i-1)/2+j-i] = sqrt(sum(WS[cbind(VV[i,],VV[j,],1:TT)]^2))
    }
  }
  class(DD) = 'dist'
  attr(DD,'Size') = NN
  return(DD)
}
