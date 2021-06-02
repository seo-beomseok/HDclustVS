#' Clustering Latent State Configurations by BFS Search.
#'
#' This function finds the cluster indices from HMM-VB or GMM-VB output.
#' @param X A data matrix.
#' @param fit HMM-VB or GMM-VB output.
#' @param C The number of clusters.
#' @return A vector of cluster indices.
#' @examples 
#' # Data generation
#' set.seed(1)
#' dat = genData2(n=300,p1=100,p2=100,C=5,rep=1) 
#' X = dat[[1]]$X_total
#' Y = dat[[1]]$z
#' n = dim(X)[1];p = dim(X)[2]
#' # The number of clusters.
#' C = 5
#' 
#' # Maximum block size is set 5% of the total dimension.
#' max.vb.size = 10
#' # Calculate Mutual Information.
#' pwmi = pairwiseMI(X)
#' # Variable block construction by ESS-DFS.
#' vbs = constVB(X,pwmi,max.vb.size)
#' 
#' # Fitting HMM-VB with 5 components.
#' fit = fitHmmvb(X,C,vbs)
#' 
#' # Semi-clusters
#' semi.cls = semicls(X,fit)
#' # Variable block selection by a bimodality test.
#' chosen.vb = semi.cls$chosen.vb
#' # Reduce the model structure
#' red.dat = reduceVB(X,fit,chosen.vb)
#' red.X = red.dat$X
#' red.vbs = red.dat$vbs
#' 
#' # Re-estimation of the HMM-VB model with reduced dimensions
#' re.fit = fitHmmvb(red.X,C,red.vbs)
#' # Final clustering
#' final.cls = finalcls(red.dat$X,re.fit,C,5)
#' @export

########################################################
####  final clusters
########################################################

finalcls <- function(X,fit,C,min.size.cls=5){
  
  DD = fit$DD
  
  tree = hclust(DD,method="complete")
  cls.ix = search.cls(tree,DD,C,min.size.cls)
  
  return(cls.ix)
}
