#' Hidden Markov Model on Variable Blocks (HMM-VB).
#'
#' This function fitt the model HMM-VB on a given data set.
#' @param X A data matrix.
#' @param C The number of clusters.
#' @param vbs Variable block structure. The output of constVB.
#' @return A list of estimated parameters for HMM-VB.
#' @import HDclust
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
####  fit HMM-VB
########################################################

fitHmmvb <- function(X,C,vbs,numst=rep(0,length(vbs$vb))){
  #-------------------------
  dim = dim(X)[2]
  nb = length(vbs$vb)
  bdim = table(vbs$ix)
  numst = ifelse(numst==0,rep(C,nb),numst)
  varorder = vbs$vb
  #-------------------------
  VbStructure <- vb(nb, dim, bdim, numst, varorder)
  hmmvb <- hmmvbTrain(X, VbStructure, trControl=trainControl(diagCov=T))
  clust <- hmmvbClust(X, model=hmmvb,
                      control = clustControl(minSize=1, modeTh = 0.0, useL1norm = T))
  
  DD = distDB2(X,hmmvb,clust)
  
  KK = hmmvb@VbStructure@numst
  TT = hmmvb@VbStructure@nb
  NN = dim(X)[1]
  BD = hmmvb@VbStructure@bdim     # size of each block
  HC = hmmvb@HmmChain             # hmm chain
  VV = t(matrix(unlist(clust@clustParam$vseq[clust@clsid]),TT))+1
  
  return(list(DD=DD,KK=KK,TT=TT,NN=NN,HC=HC,VV=VV,BD=BD,C=C,vbs=vbs,out=list(hmmvb=hmmvb,clust=clust)))
}