#' Semi-clusters
#'
#' This fuction computes semi-clusters.
#' @param X A data matrix.
#' @param fit HMM-VB or GMM-VB output.
#' @param alpha A false discovery rate (FDR) in bimodality tests of normalized mutual information (NMI) between semi-clusters and latent space states. Higher alpha level increases the chance to reduce variable blocks, but how many variable blocks will be removed depends on the bimodal distribution of NMI's and is hardly affected by alpha.
#' @return A list of semi-clusters and chosen variable block index.
#' \item{semi.cls.ix}{Semi-cluster indices.}
#' \item{chosen.vb}{Chosen relevant variable block index set.}
#' \item{normalized.mi}{Normalized mutual information between the semi-clusters and the latent states of variable blocks.}
#' \item{min.size}{Minimum size of formed semi-clusters (eta).}
#' \item{num.of.cls}{The number of formed semi-clusters (gamma).}
#' \item{pvalue}{FDR corrected p-value of the bimodality test of NMI between latent states and the chosen semi-clusters}
#' @importFrom mclust Mclust
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
####  compute semi-clusters
########################################################

semicls <- function(X,fit,alpha=0.3, lower.mode.threshold=1, verbose=FALSE){
  
  KK = fit$KK
  TT = fit$TT
  NN = fit$NN
  
  DD = fit$DD
  
  HC = fit$HC
  VV = fit$VV
  BD = fit$BD
  
  unique.KK = apply(VV,2,function(x) length(unique(x)))  
  
  tree = hclust(DD,method="complete")
  eta = 2
  CC = NN
  progress_bar = c('-','\\','|','/')
  semi.cls.list = list()
  pvalue.list = list()
  mean1.list = list()
  it = 1
  pre.cls.ix = rep(0,dim(X)[1])

  while(1){
    cls.ix = search.cls(tree,DD,NULL,eta)
    
    CC = length(table(cls.ix))
    
    if(CC<=2){
      chosen.vb = 1:TT
      break
    }
    
    if(all(cls.ix==pre.cls.ix)){
      eta = eta+1
      next
    }else{
      pre.cls.ix = cls.ix
    }

    mi.cls = rep(NA,TT)
    je.cls = rep(NA,TT)
    for(tt in 1:TT){
      if(unique.KK[tt]==1){
        mi.cls[tt] = 0
        je.cls[tt] = 1
        normalized.mi = mi.cls/je.cls
      } else{
        y2d = discretize2d(as.numeric(factor(cls.ix)),
                                    as.numeric(factor(VV[,tt])),
                                    numBins1=CC, numBins2=unique.KK[tt])
        mi.cls[tt] = mi.empirical(y2d)
        je.cls[tt] = entropy(y2d)
        normalized.mi = mi.cls/je.cls
      }
    }
    
    gmm1 = Mclust(normalized.mi,G=1,verbose=FALSE)
    gmm2 = Mclust(normalized.mi,G=2,verbose=FALSE)
    pvalue = 1-pchisq(gmm2$loglik-gmm1$loglik, df=gmm2$df-gmm1$df) 
    
    decision = rep(1,TT)
    decision[which(gmm2$classification==which.max(gmm2$parameters$mean) & normalized.mi>min(gmm2$parameters$mean))] = 2
    
    mean1 = mean(normalized.mi[decision==1])
    
    if(verbose){
      print(paste('eta=',eta,' / CC=',CC,' / pvalue=',round(pvalue,10),' / mean1=',round(mean1,3),sep=''))
      hist(normalized.mi,30)
    }
    
    chosen.vb = which(decision==2)
    
    semi.cls.list[[it]] = list(semi.cls.ix=cls.ix, 
                                  chosen.vb = chosen.vb, 
                                  normalized.mi=normalized.mi,
                                  min.size = eta,
                                  num.of.cls = CC,
                                  pvalue = pvalue,
                                  mean1 = mean1)
    pvalue.list[[it]] = pvalue
    mean1.list[[it]] = mean1
    it = it+1
    
    eta = eta+1
    cat(paste(progress_bar[eta%%4+1],'\b',sep=""))
    
  }
  
  qvalue.list = p.adjust(pvalue.list,method='BH')
  if(any((qvalue.list<=alpha) & (mean1.list<=lower.mode.threshold))){
    j1 = which(mean1.list<=lower.mode.threshold)
    j2 = which.min(qvalue.list[j1])
    j = j1[j2]
    semi.cls = semi.cls.list[[j]]
  }else{
    j = which.min(qvalue.list)
    semi.cls = semi.cls.list[[j]]
    semi.cls$chosen.vb = 1:TT
  }
  
  return(semi.cls)
}
#---------------------------------
