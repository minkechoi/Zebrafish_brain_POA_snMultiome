#' creates an `eulerr::`-friendly vector to plot overlaps of two sets.
calc_overlaps <- function(x){
  a = length(which(!(x[[1]] %in% x[[2]])))
  b = length(which(!(x[[2]] %in% x[[1]])))
  a_and_b = length(which(x[[2]] %in% x[[1]]))
  
  v = c(a,b,a_and_b)
  n = names(x)
  
  res = setNames(v, c(n, paste(n, collapse = "&") ) )
  
  return(res)
}
