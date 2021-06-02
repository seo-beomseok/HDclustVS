#' Computing the size of the subtree of a node.
#'
#' This function calculates the size of leaves-set of a node.
#' @param M A matrix of merges of nodes, which is obtained from hclust() output.
#' @return A vector with the size of leaves-set of each node.
#' @export

########################################################
####  node size
########################################################

nodeSize = function(M){
  p = dim(M)[1]+1
  SZ = rep(NA,p-1)
  for(i in 1:(p-1)) {
    if(M[i,1]<0 && M[i,2]<0) SZ[i] = 2
    else if(M[i,1]<0 && M[i,2]>0) SZ[i] = SZ[M[i,2]]+1
    else if(M[i,1]>0 && M[i,2]<0) SZ[i] = SZ[M[i,1]]+1
    else if(M[i,1]>0 && M[i,2]>0) SZ[i] = SZ[M[i,1]]+SZ[M[i,2]]
  }
  return(SZ)
}