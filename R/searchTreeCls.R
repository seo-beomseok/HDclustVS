#' Clustering Latent State Configurations by BFS Search.
#'
#' This function searches clusters by BFS search from a dendrogram tree.
#' @param tree An object of hclust result.
#' @param dist The distance measures between nodes for the tree.
#' @param num.cls The number of cluster. If it is NULL, it finds semi-cluster. If it is specified as an integer, it finds clusters with the number of groups.
#' @param min.size The minimum size of each cluster. Default is 5. 
#' @return A vector of cluster indices for each leaf node.
#' DD = dist(X)
#' min.size.cls = 5
#' tree = hclust(DD)
#' cls.ix = search.cls(tree,DD,C,min.size.cls)
#' @export

########################################################
####  search clusters from a tree
########################################################

search.cls = function(tree,dist,num.cls=NULL,min.size=5){
  ind = NULL
  cand = NULL
  
  NN = nrow(tree$merge)+1
  SZ = nodeSize(tree$merge)
  
  search.cand.node = function(i){
    m = tree$merge[(NN-1):(NN-i),]
    
    cand = m[which(m<(NN-i))]
    id = c(SZ,1)[ifelse(cand<0,NN,cand)]>=min.size
    ind = cand[id]
    return(ind)
  }
  
  num.cl = rep(NA,NN-1)
  for(i in 1:NN){
    ind = search.cand.node(i)
    num.cl[i] = length(ind)
  }
  
  if(is.null(num.cls)){
    ind = search.cand.node(which.max(num.cl))
  }else if(length(which(num.cl==num.cls))==0){
    return(NA)
  }else{
    ind = search.cand.node(min(which(num.cl==num.cls)))
  }
  
  cls = lapply(1:length(ind), function(i) search.leaves(tree$merge,ind[i]))
  clsd = unlist(cls)
  unclsd = (1:NN)[-clsd]
  distmat = as.matrix(dist)
  mm = lapply(cls,function(i) distmat[i,unclsd])
  
  if(length(unclsd)>1){
    unclid = apply(sapply(mm,function(i) apply(i,2,min)),1,which.min)
  }else if(length(unclsd)==1){
    unclid = which.min(sapply(mm,min))
  }else{
    unclid = NULL
  }
  
  clsix = rep(NA,NN)
  for(i in 1:length(cls)){
    clsix[cls[[i]]] = i
  }
  clsix[unclsd] = unclid
  return(clsix)
}
