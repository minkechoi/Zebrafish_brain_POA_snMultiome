#' @params x: a counts object from seurat; ident.e. mydata@assays$RNA@counts
#' @params identities: a (named) vector with cluster information on each cell. 
#' IMPORTANT: cell order is assumed the same between x and identities.
pseudobulk <- function(x, identities){
  
  if (!is.factor(identities)) identities <- factor(identities, levels = unique(identities))
  
  idents <- levels(identities)
  
  mtx <- matrix(0, nrow = nrow(x), ncol = length(idents))
  rownames(mtx) <- rownames(x)
  colnames(mtx) <- idents
  
  for (ident in idents){
    filt <- which(identities == ident)
    x_i <- x[,filt]
    pseudocounts_cluster_i <- rowSums(x_i)
    mtx[,ident] <- pseudocounts_cluster_i
    message("Done cluster ",ident)
  }
  
  return(mtx)
}