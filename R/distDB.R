#' # Davis-Bouldin Type Distance
#'
#' This function calculates Davis-Bouldin type distance.
#' @param mean A matrix of mean values.
#' @param Sigma A list of variance-covariance matrix.
#' @return A distance matrix.
#' @export

########################################################
####  Davis-Bouldin type distance
########################################################

distDB = function(mean,Sigma){
  WW = mean
  SS = Sigma
  
  NN = dim(WW)[1]
  DD = rep(NA,NN*(NN-1)/2)
  for(i in 1:(NN-1)){
    for(j in (i+1):NN){
      DD[NN*(i-1)-i*(i-1)/2+j-i] = sqrt(sum(((WW[i,]-WW[j,])^2)/(SS[i,]+SS[j,])))
    }
    #print(i)
  }
  
  class(DD) = 'dist'
  attr(DD,'Size') = NN
  return(DD)
}