#' Variable Blocks Construction by ESS-DFS.
#'
#' This function constructs variable blocks by ESS-DFS algorithm.
#' @param X A data matrix.
#' @param pwmi A matrix of pairwise mutual information.
#' @param max.size Maximum size of variable blocks.
#' @return A list of the variable indices grouped by the constructed variable blocks.
#' \item{vb}{variable indices for variable blocks.}
#' \item{ix}{variable block indices for variables.}
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
#' # Variable block construction by ESS-DFS.
#' vbs = constVB(X,pwmi,max.vb.size)
#' @export

########################################################
####  Construct Variable Blocks
########################################################


constVB <- function(X,pwmi,max.size = dim(X)[2]/10){
  
  DD = as.dist(t(1/pwmi))
  hclmi = hclust(DD, method = 'complete')
  
  M = hclmi$merge
  O = hclmi$order
  H = hclmi$height
  
  node.size = nodeSize(M)
  
  candi.node = which(node.size <= max.size)
  fst.rootnode = candi.node[which.max(search.parent(M,candi.node))]  # find initial root node 
  
  seq.rootnode = search.rootnode(M,fst.rootnode,max.size,DD)         # root node index of each block
  
  vb = search.vb(M,seq.rootnode)
  
  return(vb)
}