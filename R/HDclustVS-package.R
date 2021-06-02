#' Block-wise Variable Selection for Clustering via Latent States of Mixture Models
#' 
#' \pkg{HDclustVS} is an R package for new block-wise variable selection methods for clustering, HMM-VB-VS and GMM-VB-VS, which exploit the latent states of the hidden Markov model on variable blocks or Gaussian mixture model. The variable blocks are formed by early-stop-and-sorted-depth-first-search (ESS-DFS) on a dendrogram created based on the mutual information between any pair of variables. Then, the variable selection is conducted by an independence test between the latent states and semi-clusters which are the smaller clusters that will be further grouped into final clusters.
#' @aliases HDclustVS
#' @aliases HDclustVS-package
#' @keywords internal
"_PACKAGE"
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