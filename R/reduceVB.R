#' Dimension-reduced Data Set.
#'
#' This fuction computes the dimension-reduced data set.
#' @param X A data matrix.
#' @param fit HMM-VB or GMM-VB output.
#' @param chosen.vb A vector of chosen variable block indices.
#' @return A list of dimension-reduced data set and its variable block structure.
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
####  return reduced data and reduced variable blocks
########################################################

reduceVB <- function(X,fit,chosen.vb){
  A = chosen.vb
  #-----------------
  X = X
  vb = fit$vbs$vb
  ix = fit$vbs$ix
  #-----------------
  chosen.X = (1:dim(X)[2])[unlist(vb[A])]
  X.A = X[,chosen.X]
  vb.A = vb[A]
  ix.A = rep(NA,dim(X.A)[2])
  for(i in 1:length(A)){
    vb.A[[i]] = match(vb.A[[i]],chosen.X)
    ix.A[vb.A[[i]]]=i
  }
  #-----------------
  return(list(X=X.A,vbs=list(vb=vb.A,ix=ix.A)))
}
