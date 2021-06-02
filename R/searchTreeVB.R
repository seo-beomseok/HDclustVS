#' Search Tree Graph for Variable Block Construction.
#'
#' Functions used to construct variable blocks.
#' @name searchTreeVB
#' @aliases search.parent
#' @aliases search.leaves
#' @aliases search.rootnode
#' @aliases search.vb
#' @param M A matrix of merges which is obtained from hclust() output.
#' @param loc An index or a sequence of indices of nodes in the tree.
#' @param max.size Maximum size of variable blocks.
#' @param dist.measure Distance measure of the tree.
#' @param seq A sequence of root nodes of subtrees corresponding to variable blocks
#' @examples 
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
#' 
#' DD = as.dist(t(1/pwmi))
#' hclmi = hclust(DD, method = 'complete')
#' M = hclmi$merge
#' node.size = nodeSize(M)
#' candi.node = which(node.size <= max.vb.size)
#' fst.rootnode = candi.node[which.max(search.parent(M,candi.node))]  # find initial root node 
#' seq.rootnode = search.rootnode(M,fst.rootnode,max.vb.size,DD)         # root node index of each block
#' vb = search.vb(M,seq.rootnode)

########################################################
####  search parent node
########################################################
#' @rdname searchTreeVB
#' @return search.parent gives a vector of parent node indices of given node indices in a tree.
#' @export

search.parent = function(M,loc){
  p = nrow(M)+1
  if(min(loc) < 1 | max(loc) >= p) stop("loc must be between 1 and n-1")
  return(match(loc,M)%%nrow(M))
}

########################################################
####  search leaves node
########################################################
#' @rdname searchTreeVB
#' @return search.leaves gives a vector of leaves sets of a given node index in a tree.
#' @export

search.leaves = function(M,loc){
  p = nrow(M)+1
  
  if(loc < 0 & loc >= -p) return(-loc)
  else if(loc < -p | loc >= p) stop("loc must be between -n and n-1")
  
  ind = rep(NA,p) # individual num list
  cc = rep(NA,p)  # chosen next path
  dd = rep(NA,p)  # unchosen path
  
  i=1 # index to save in ind[]
  j=1 # index for present path in cc[]
  l=1 # index to save in cc[]
  m=1 # index to save in dd[]
  
  cc[1]=loc
  
  while(1){
    while(!is.na(cc[j])){
      if(M[cc[j],1]>0 && M[cc[j],2]>0) {
        cc[l+1] = max(M[cc[j],1],M[cc[j],2])
        dd[m+1] = min(M[cc[j],1],M[cc[j],2])
        l=l+1
        m=m+1
      } else if(M[cc[j],1]<0 && M[cc[j],2]>0) {
        ind[i] = -M[cc[j],1]
        i=i+1
        cc[l+1] = M[cc[j],2]
        l=l+1
      } else if(M[cc[j],1]>0 && M[cc[j],2]<0) {
        ind[i] = -M[cc[j],2]
        i=i+1
        cc[l+1] = M[cc[j],1]
        l=l+1
      } else if(M[cc[j],1]<0 && M[cc[j],2]<0) {
        ind[i] = max(-M[cc[j],1],-M[cc[j],2])
        ind[i+1] = min(-M[cc[j],1],-M[cc[j],2])
        i=i+2
      }
      j=j+1
    }
    cc = rep(NA,p)
    j=1
    l=1
    if(sum(which(!is.na(dd)))>0) {
      cc[1] = dd[max(which(!is.na(dd)))]
      dd[max(which(!is.na(dd)))] = NA
    } else break
  }
  
  return(ind[1:max(which(!is.na(ind)))])
}

########################################################
####  search root node of each variable blocks
########################################################
#' @rdname searchTreeVB
#' @return search.rootnode gives a vector of root node indices for variable blocks given the first node (loc), maximum size of block (max.size), and distance measure (dist).
#' @export

search.rootnode <- function(M,loc,max.size,dist.measure,subgraph=F) {
  
  DD = dist.measure
  SZ = nodeSize(M)
  
  p = nrow(M)+1
  ind = rep(NA,p) # individual num list
  cc = rep(NA,p)  # chosen next path
  dd = rep(NA,p)  # unchosen path
  ee = NA # highest node which have been examined
  
  i=1 # index to save in ind[]
  j=1 # index for present path in cc[]
  l=1 # index to save in cc[]
  m=1 # index to save in dd[]
  
  cc[1]=ee=loc
  vloc = which(M==ee)
  dd[1]=ifelse(vloc>(p-1),M[vloc-(p-1),1],M[vloc,2])
  
  while(1){
    
    while(!is.na(cc[j])){
      if(cc[j]<0) {
        ind[i] = cc[j]
        i=i+1
      } else if(SZ[cc[j]]<=max.size) { 
        ind[i] = cc[j]
        i=i+1
      } else if(SZ[cc[j]]>max.size) {
        last.ind = ind[max(which(!is.na(ind)))]
        s0 = search.leaves(M,last.ind)
        s1 = search.leaves(M,M[cc[j],1])
        s2 = search.leaves(M,M[cc[j],2])
        
        dmin = which.min(c(max(as.matrix(DD)[s0,s1]),max(as.matrix(DD)[s0,s2])))
        dmax = 3-dmin
        
        cc[l+1] = M[cc[j],dmin]
        dd[m+1] = M[cc[j],dmax]
        l=l+1
        m=m+1
      } 
      j=j+1
    }
    
    if(subgraph==T) dd[1] = NA
    dd.ix = which(!is.na(dd))
    
    if(sum(dd.ix)==0) {
      if(subgraph==T){
        return(ind[1:max(which(!is.na(ind)))])
        break
      }
      vloc0 = which(M==ee)  
      vloc1 = which(M==ifelse(vloc0>(p-1),vloc0-(p-1),vloc0))
      if(sum(vloc1)!=0) {
        cc = rep(NA,p)
        cc[1] = ee = ifelse(vloc1>(p-1),M[vloc1-(p-1),1],M[vloc1,2])
        j=1
        l=1
      } else break
    } else if(sum(dd.ix)>0) {
      cc = rep(NA,p)
      j=1
      l=1
      cc[1] = dd[m]
      dd[m] = NA
      m=m-1
    }
    
  }
  
  return(ind[1:max(which(!is.na(ind)))])
}


########################################################
####  search leaves of each variable blocks
########################################################
#' @rdname searchTreeVB
#' @return search.vb gives a list of leaves sets for each root node (seq).
#' @export

search.vb = function(M,seq){
  p = nrow(M)+1
  
  cls = list()
  clsix = rep(NA,p)
  
  for(m in 1:length(seq)) {
    if(seq[m]>0){
      cls[[m]] = search.leaves(M,seq[m])
      clsix[cls[[m]]] = m
    } else if(seq[m]<0){
      cls[[m]] = -seq[m]
      clsix[cls[[m]]] = m
    }
  }
  return(list(vb=cls,ix=clsix))
}

