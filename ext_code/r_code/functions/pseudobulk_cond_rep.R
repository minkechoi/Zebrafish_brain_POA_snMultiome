#' @params x: a counts object from seurat; ident.e. mydata@assays$RNA@counts
#' @params identities: a (named) vector with cluster information on each cell. 
#' @params conditions: a (named) vector with experiment information on each cell. 
#' @params replicates: a (named) vector with replicate information on each cell. 
#' IMPORTANT: cell order is assumed the same between x and identities, conditions, and replicates.
pseudobulk_cond_rep <- function(x, identities, conditions, replicates){
  
  # parse the different filtering vectors: turn into factors
  if (!is.factor(identities)) identities <- factor(identities, levels = unique(identities))
  if (!is.factor(conditions)) conditions <- factor(conditions, levels = unique(conditions))
  if (!is.factor(replicates)) replicates <- factor(replicates, levels = unique(replicates))
  
  # retrieve unique instances
  idents <- levels(identities); message("found idents")
  conds <- levels(conditions); message("found conds")
  reps <- levels(replicates); message("found reps")
  
  # create matrix where to put counts
  mtx <- 
    matrix(
      0, 
      nrow = nrow(x), # all same genes as in input data
      ncol = length(idents)*length(conds)*length(reps) # amount of combinations of cluster, experiment, replicate
    )
  rownames(mtx) <- rownames(x) # all same genes as in input data
  colnames(mtx) <- apply(expand.grid(idents,conds, reps), 1, paste, collapse="_") # combinations of cluster, experiment, replicate
  
  # create a data frame with information for DGE
  sampletable <- data.frame(
    sample = paste0("S",formatC(1:ncol(mtx),width=2,flag="0")),
    id_combined = apply(expand.grid(idents,conds, reps), 1, paste, collapse="_"),
    ctype = "",
    condition = "",
    replicate = "",
    ncells = 0
  )
  
  # begin nested loop
  for (ident in idents) {
    for (cond in conds) {
      for (rep in reps) {
        sample <- paste(ident,cond,rep,sep="_") # recreate a sample name
        message("Starting sample ",sample)
        
        filt = identities == ident & conditions == cond & replicates == rep # combination of conditions: cells from cluster i in experiment j replicate z
        
        ncells_filt = length(which(filt)==TRUE)
        
        message("Found ",ncells_filt," cells of sample ",sample)
        
        if(ncells_filt == 0) {
          message("No cells in ", sample, ", skipping")
          next 
        } else {
          x_i <- x[,filt] # filter the matrix of counts x cells, keep cells from combination filter
        }
        
        if(ncells_filt > 1){ # if we have more than one cell passing the filter
          pseudocounts_sample_i <- rowSums(x_i) # sum the counts of all the same cells, output at a gene level
        } else {
          pseudocounts_sample_i <- x_i # otherwise this column will be the columns of that single cell
        }
        
        which_sample <- which(colnames(mtx) == sample)
        mtx[,which_sample] <- pseudocounts_sample_i #  dump those values in the corresponding sample column of the output matrix
        sampletable[which_sample,3:5] <- c(ident,cond,rep) # dump the values of sample etc in the corresponding row of sampletable data
        sampletable[which_sample,6] <- ncells_filt
        message("Done sample ",sample) #  message to screen
      }
    }
  }
  
  colnames(mtx) <- sampletable$sample[match(sampletable$id_combined,colnames(mtx))]
  
  # create results list
  res <- list(
    matrix = mtx, # output matrix
    sampletable = sampletable # output data frame with information for DGE
  )
  
  return(res) # return the results list
}