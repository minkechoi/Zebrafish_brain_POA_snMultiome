deseq_pseudobulk <- 
  function(
    count_matrix, samples_info, celltype, ...
  ) {
    
    if( !(celltype %in% samples_info$ctype)) 
      stop("E: cell type specified is not in matrix count or samples table.")

    # Parse and filter input data
    cell <- celltype
    cell_filt = which(samples_info$ctype == cell)

    # filter count matrix
    m <- count_matrix[,cell_filt]
    rownames(m) <- rownames(count_matrix)
    m <- m[rowSums(m) > 1,]
    
    # filter samples info
    d <- samples_info[cell_filt,]
    
    res <- deseq_sc(m = m, d = d, cell = cell, ...)
    
    return(res)
  }
