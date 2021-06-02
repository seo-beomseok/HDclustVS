#' Cross Validation for Variable Selection Threshold
#'
#' This function finds a threshold for variable selection using cross validation.
#' @name viterbi
#' @aliases TransProb
#' @aliases EmissProb
#' @aliases Viterbi
#' @aliases GaussParam
#' @param X A data matrix
#' @param fit HMM-VB or GMM-VB output
#' @param TP Transition probability list. The output of TransProb().
#' @param EP Emission probability list. The output of EmissProb().
#' @param NN Sample size. dim(X)[1].
#' @param TT The length of variable blocks.
#' @param KK A vector of the number of components for variable blocks.
#' @param HC Hidden Markov chain
#' @param VV Viterbi path. The output of Viterbi().
#' @param BD A vector of the number of variables for each variable block.

########################################################
####  Cross Validation for Variable Selection Threshold
########################################################
#' @rdname viterbi
#' @return TransProb gives a list of transition probability.
#' @export

TransProb = function(X,fit){
  
  TT = fit$TT                     # number of blocks
  NN = dim(X)[1]                  # number of observations
  
  HC = fit$HC                     # hmmChain
  
  TP = list()               # Transition probability
  for(t in 1:TT) TP[[t]] = HC[[t]]@a
  
  return(TP)
}

#' @rdname viterbi
#' @return EmissProb gives a list of emission probability.
#' @export

EmissProb = function(X,fit){
  
  KK = fit$KK                     # number of components for each block
  TT = fit$TT                     # number of blocks
  NN = dim(X)[1]                  # number of observations
  
  HC = fit$HC                     # hmmChain
  VB = fit$vbs$vb                 # variable blocks
  
  EP = list()               # Emission probability
  for(t in 1:TT){
    EP[[t]] = list()
    
    for(n in 1:NN){
      EP[[t]][[n]] = rep(NA,KK[[t]])
      
      ME = HC[[t]]@mean
      SIG = HC[[t]]@sigma
      for(k in 1:KK[[t]]){
        EP[[t]][[n]][k] = mvtnorm::dmvnorm(X[n,VB[[t]]],mean=ME[k,],sig=SIG[[k]])
      }
      EP[[t]][[n]] = EP[[t]][[n]]/max(EP[[t]][[n]])
    }
  }
  
  {                    # When EP generates NaN
    for(n in 1:NN){
      for(t in 1:TT){
        for(k in 1:KK[t])
          if(is.nan(EP[[t]][[n]][k])){
            EP[[t]][[n]][k]=0
            warning("Emmision Probability generates NaN")
          }
        }
      }
  }
  
  return(EP)
}
  
#' @rdname viterbi
#' @return Viterbi gives a matrix of viterbi paths.
#' @export
  
Viterbi <- function(TP,EP,NN,TT,KK){
  
  #llik = rep(NA,NN)
  VV = matrix(NA,NN,TT)     # Viterbi path
  for(n in 1:NN){
    
    T1 = matrix(NA,max(KK),TT)
    T2 = matrix(NA,max(KK),TT)
    
    T1[,1] = TP[[1]]*EP[[1]][[n]]
    T2[,1] = 0
    
    for(t in 2:TT){
      T1[,t] = apply(t(T1[,t-1]*TP[[t]])*EP[[t]][[n]],1,max)
      T2[,t] = apply(t(T1[,t-1]*TP[[t]])*EP[[t]][[n]],1,which.max)
    }
    
    x = rep(NA,TT)
    
    x[TT] = which.max(T1[,TT])
    for(t in TT:2){
      x[t-1] = T2[x[t],t]
    }
    VV[n,] = x
    
    #llik[n] = sum(log(T1[cbind(x,1:TT)]))
  }
  return(VV)
}
  
#' @rdname viterbi
#' @return GaussParm gives a list of estimated mixture model parameters for each sample.
#' \item{WW}{Mean values of each sample corresponding to MAP state configuration.}
#' \item{SS}{Variance values of each sample corresponding to MAP state configuration.}
#' @export

GaussParam <- function(HC,VV,BD,NN,TT,KK){
  WW = matrix(NA,NN,sum(BD))
  SS = matrix(NA,NN,sum(BD))
  for(n in 1:NN){
    for(t in 1:TT){
      WW[n,(sum(BD[1:t-1])+1):sum(BD[1:t])]=HC[[t]]@mean[VV[n,t],]
      SS[n,(sum(BD[1:t-1])+1):sum(BD[1:t])]=diag(HC[[t]]@sigma[[VV[n,t]]])
    }
  }
  return(list(WW=WW,SS=SS))
}


  