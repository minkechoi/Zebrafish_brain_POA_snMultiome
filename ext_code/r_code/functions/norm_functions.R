# Define custom normalization methods

#' Function adapted from Booeshaghi et al., 2021
do_pf <- function(mtx, sf = NULL) {
  require(Matrix)
  pf <- Matrix(rowSums(mtx), ncol = 1) # sum of counts per cell
  if (is.null(sf)) {
    sf <- mean(pf) # calculate scaling factor
  }
  dsf <- Diagonal(x = as.numeric(sf/pf)) # diagonal sparse matrix with sf/pf values on the diagonal
  res <- dsf %*% mtx # matrix multiplication should return same dimensions as cts
  rownames(res) <- rownames(mtx)
  colnames(res) <- colnames(mtx)
  
  return(res)
}

PFlog1pPF <- function(mtx) {
  return(
    do_pf(Matrix(log1p(do_pf(mtx))))
  )
}

#' Function adapted from Booeshaghi et al., 2021, the scaling factor is calculated from a specified quantile q
do_pf_q <- function(mtx, q = 0.5, sf = NULL) {
  pf <- Matrix(rowSums(mtx), ncol = 1) # sum of counts per cell
  if (is.null(sf)) {
    sf <- quantile(pf,q) # calculate scaling factor
  }
  dsf <- Diagonal(x = as.numeric(sf/pf)) # diagonal sparse matrix with sf/pf values on the diagonal
  res <- dsf %*% mtx # matrix multiplication should return same dimensions as cts
  rownames(res) <- rownames(mtx)
  colnames(res) <- colnames(mtx)
  
  return(res)
}

#' "Cell weight" Matrix (we need a better name). Developed in collaboration with Isabel Liao from Yi-Jyun's Lab
#' @params x: a gene per cell cluster matrix with number of counts per gene
#' @params y: a gene per cell cluster matrix with how many cells are expressing a given gene per cluster
#' @params C: a named vector with the size of all clusters
#' @params min_counts: an integer with the minimum counts of a gene in the dataset
#' @params min_cells: an integer with the minimum cells expressing a gene in the dataset
get_cellweight_matrix <- function(x,y,C,min_counts = 1,min_cells=1){
  
  x = x[rowSums(x) >= min_counts,]
  
  y = y[rownames(y) %in% rownames(x),]
  y[y < min_cells] <- 0
  y = y[rowSums(y)>0,]
  
  m <- 
    matrix(
      0, 
      nrow = nrow(y),
      ncol = ncol(y),
      dimnames=dimnames(y)
    )
  
  for (j in 1:ncol(m)){
    
    a = y[,j]/C[j]
    
    b = rowSums(y[,-j])/sum(C[-j])
    
    m[,j] = a/b
  }
  
  m = 1 - exp(-m)
  
  m[is.nan(m)] <- 0
  
  return(m)
}
