#' Visulization of a partition on 2 dimensional space with convex hull.
#'
#' This function plots a partition on 2 dimensional principal components space with convex hull.
#' @param X A matrix of data.
#' @param z Cluster labels of data.
#' @return None
#' @importFrom dplyr '%>%'
#' @import ggplot2
#' @examples
#' # Data generation
#' set.seed(1)
#' dat = genData2(n=300,p1=100,p2=100,C=5,rep=1) 
#' X = dat[[1]]$X_total
#' Y = dat[[1]]$z
#' n = dim(X)[1];p = dim(X)[2]
#' 
#' PCA = prcomp(dat[[1]]$X_relev)
#' Xr = prcomp(dat[[1]]$X_relev)$x[,1:2]
#' plotCls(Xr,Y,title="Relevant var.")
#' Xi = prcomp(dat[[1]]$X_irrel)$x[,1:2]
#' plotCls(Xi,Y,title="Irrelevant var.")
#' Xt = prcomp(dat[[1]]$X_total)$x[,1:2]
#' plotCls(Xt,Y,title="Total var.")
#' @export

########################################################
####  plot with convex hull
########################################################

plotCls <- function(X,z,title="",xlab="",ylab="",legend.title="",no.legend=F){
  
  if(dim(X)[2]==2){
    Y = data.frame(PC1=X[,1],PC2=X[,2])
  } else if(dim(X)[2]>2){
    Y = as.data.frame(prcomp(X)$x)[,1:2]
  }
  
  A = as.character(z)
  
  s = cbind(Y,A) %>% split(z)
  ch = s %>% lapply(., function(el) chull(el$PC1, el$PC2))
  ch = lapply(names(ch), function(x) s[[x]][ch[[x]],]) %>% do.call(rbind, .) 
  
  mu_x = sapply(1:max(z),function(i){mean(Y[which(z==i),1])})
  mu_y = sapply(1:max(z),function(i){mean(Y[which(z==i),2])})
  
  if(no.legend){
    ggplot2::ggplot(as.data.frame(Y[,1:2]), ggplot2::aes(x=PC1, y=PC2, color = A)) +
      ggplot2::geom_point(shape=19, size=1, alpha=0.5) +
      ggplot2::geom_polygon(data = ch, ggplot2::aes(fill = A), alpha = 0.2) +
      ggplot2::theme_classic() +
      ggplot2::labs(x=xlab,y=ylab) +
      ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(face="bold", hjust = 0.5), legend.position = "none", text = element_text(size=20), plot.margin=margin(10,5,5,10)) +
      #ggplot2::guides(fill=ggplot2::guide_legend(title=legend.title), color=ggplot2::guide_legend(title=legend.title)) +
      ggplot2::annotate("text", x=mu_x, y=mu_y, label=1:max(z), size=6, fontface = "bold")
  }else{
    ggplot2::ggplot(as.data.frame(Y[,1:2]), ggplot2::aes(x=PC1, y=PC2, color = A)) +
      ggplot2::geom_point(shape=19, size=1, alpha=0.5) +
      ggplot2::geom_polygon(data = ch, ggplot2::aes(fill = A), alpha = 0.2) +
      ggplot2::theme_classic() +
      ggplot2::labs(x=xlab,y=ylab) +
      ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(face="bold", hjust = 0.5), text = element_text(size=20), plot.margin=margin(10,5,5,10)) +
      ggplot2::guides(fill=ggplot2::guide_legend(title=legend.title), color=ggplot2::guide_legend(title=legend.title)) +
      ggplot2::annotate("text", x=mu_x, y=mu_y, label=1:max(z), size=6, fontface = "bold")
  }
}
