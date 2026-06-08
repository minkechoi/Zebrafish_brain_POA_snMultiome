#' Staircase TF expr
highest_val <- function(x){ # x: a vector of named values
  a = x - mean(x)
  b = which(a == max(a) )
  if(length(b) > 1) b = b[1]
  c = names(x)[b]
  return(c)
}

highest_val_0 <- function(x){ # x: a vector of named values
  a = x
  b = which(a == max(a) )
  if(length(b) > 1) b = b[1]
  c = names(x)[b]
  return(c)
}