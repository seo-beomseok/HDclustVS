#' Pairwise Mutual Information.
#'
#' This function calculates pairwise mutual information between variables.
#' @param X A data matrix.
#' @return A matrix of pairwise mutual information.
#' @importFrom entropy discretize2d
#' @importFrom entropy mi.empirical
#' @importFrom entropy entropy
#' @examples 
#' # Data generation
#' set.seed(1)
#' dat = genData2(n=300,p1=100,p2=100,C=5,rep=1) 
#' X = dat[[1]]$X_total
#' Y = dat[[1]]$z
#' n = dim(X)[1];p = dim(X)[2]
#' # The number of clusters.
#' C = 5
#' # Maximum block size is set 5% of the total dimension.
#' max.vb.size = 10
#' # Calculate Mutual Information.
#' pwmi = pairwiseMI(X)
#' @export

########################################################
####  pairwise mutual information
########################################################

pairwiseMI <-function(X,numBins1=10,numBins2=10){
  n = dim(X)[1]
  p = dim(X)[2]
  
  MI = matrix(NA,p,p)
  #JE = matrix(NA,p,p)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      y2d = discretize2d(X[,i], X[,j], numBins1, numBins2)
      MI[i,j] = mi.empirical(y2d)
      #JE[i,j] = entropy(y2d)
    }
    #print(i)
  }
  return(MI)
}


