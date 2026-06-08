#' Returns a named (=cell) list of scores.
#' Should be downstream compatible with boxplot(var~var, dataset)
gene_score <- function(x, gene_set, gene_pool = rownames(x), remove_set_from_pool = FALSE, fraction = 1, seed = 4343){
  # parse gene pool
  ## remove gene set if asked
  if (remove_set_from_pool == TRUE){
    filt <- !(gene_pool %in% gene_set)
    gene_pool <- gene_pool[filt]
  }
  
  ## subsample
  n_pool <- length(gene_pool)
  n <- trunc(n_pool*fraction) # if fraction == 1 this will be just all of them
  set.seed(seed) # user can control the seed
  pool <- sample(gene_pool, n)
  
  # calculate score
  a_i <- apply(x[gene_set,], 2, mean) # is there a faster way?
  b_i <- apply(x[pool,], 2, mean) # is there a faster way?
  score <- a_i - b_i
  
  # result
  return(score)
  
}