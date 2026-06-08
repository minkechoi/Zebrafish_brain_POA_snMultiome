#' @params x: a counts object from seurat; ident.e. mydata@assays$RNA@counts
#' @params identities: a (named) vector with cluster information on each cell.
#' @params min_counts: numeric, report cells with more than this number of reads for a given gene.
#' IMPORTANT: cell order is assumed the same between x and identities.
pseudobulk_ncells <- function(x, identities, min_counts = 1){
  
  if (!is.factor(identities)) identities <- factor(identities, levels = unique(identities))
  
  idents <- levels(identities)
  
  mtx <- matrix(0, nrow = nrow(x), ncol = length(idents))
  rownames(mtx) <- rownames(x)
  colnames(mtx) <- idents
  
  for (ident in idents){
    filt <- which(identities == ident)
    x_i <- x[,filt]
    ncells_cluster_i <- 
      apply(
        x_i,
        MARGIN = 1, 
        FUN = function(x){
          length( which( x >= min_counts ) )
        }
      )
    mtx[,ident] <- ncells_cluster_i
    message("Done cluster ",ident)
  }
  
  return(mtx)
}
