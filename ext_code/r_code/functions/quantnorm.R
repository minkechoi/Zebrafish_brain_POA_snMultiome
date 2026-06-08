#' 
#' quantnorm: prepare a quantile normalisation of input data
#' quantnorm(a) --> a_q
#' 
quantnorm <- function(x){
  require(preprocessCore)
  colnames_x <- colnames(x)
  rownames_x <- rownames(x)
  x <- as.matrix(x)
  x_q <- normalize.quantiles(x) # figure out why we get zeroes and neg numbers here
  colnames(x_q) <- colnames_x
  rownames(x_q) <- rownames_x
  return(x_q)
}