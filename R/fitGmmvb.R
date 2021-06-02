#' Gaussian Mixture Model on Variable Blocks (GMM-VB).
#'
#' This function fits the model GMM-VB on a given data set.
#' @param X A data matrix.
#' @param C The number of clusters.
#' @param vbs Variable block structure. The output of constVB.
#' @return A list of estimated parameters for GMM-VB.
#' @importFrom mclust Mclust
#' @importFrom mvtnorm dmvnorm
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
#' fit = fitGmmvb(X,C,vbs)
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
#' re.fit = fitGmmvb(red.X,C,red.vbs)
#' # Final clustering
#' final.cls = finalcls(red.dat$X,re.fit,C,5)
#' @export

########################################################
####  fit GMM-VB
########################################################

fitGmmvb <- function(X,C,vbs,numst=rep(0,length(vbs$vb))){
  #-------------------------
  dim = p = dim(X)[2]
  nb = length(vbs$vb)
  bdim = table(vbs$ix)
  numst = ifelse(numst==0,rep(C,nb),numst)
  varorder = vbs$vb
  #-------------------------
  
  KK = numst
  TT = nb
  NN = dim(X)[1]
  BD = bdim
  
  MCLT = Mclust(X,C)
  
  Mmu=MSig=list()
  for(i in 1:C){
    Mmu[[i]] = MCLT$parameters$mean[,i]
    MSig[[i]] = MCLT$parameters$variance$sigma[,,i]
  }
  
  llt = array(dim=c(NN,TT,C))
  vvt = array(dim=c(NN,TT))
  for(i in 1:NN){
    for(t in 1:TT){
      for(j in 1:numst[1]){
        if(bdim[t]==1){
          llt[i,t,j] = dnorm(X[i,varorder[[t]]],Mmu[[j]][varorder[[t]]],diag(MSig[[j]])[varorder[[t]]],log=T)
        }else{
          llt[i,t,j] = mvtnorm::dmvnorm(X[i,varorder[[t]]],Mmu[[j]][varorder[[t]]],MSig[[j]][varorder[[t]],varorder[[t]]],log=T)
        }
      }
      vvt[i,t] = which.max(llt[i,t,])
    }
  }
  
  
  WW = matrix(NA,NN,p)
  SS = matrix(NA,NN,p)
  VV = matrix(NA,NN,nb)
  for(t in 1:TT){
    VV[,t] = vvt[,t]
    WW[,varorder[[t]]] = t(MCLT$parameters$mean[varorder[[t]],vvt[,t]])
    SS[,varorder[[t]]] = sapply(varorder[[t]],function(x) MCLT$parameters$variance$sigma[x,x,vvt[,t]])
  }
  
  DD = distDB(WW,SS)
  
  HC = list()
  for(t in 1:TT){
    HC[[t]] = new('HMM',
                dim=BD[t],
                numst=KK[t],
                mean=t(MCLT$parameters$mean[varorder[[t]],]),
                sigma=lapply(1:KK[t],function(x) MCLT$parameters$variance$sigma[varorder[[t]],varorder[[t]],x]))
  }
  
  return(list(DD=DD,KK=KK,TT=TT,NN=NN,HC=HC,VV=VV,BD=BD,C=C,vbs=vbs,out=MCLT))
}
