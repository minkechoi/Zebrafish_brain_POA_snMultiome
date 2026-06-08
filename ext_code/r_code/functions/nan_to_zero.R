nan_to_zero <- function(x){
  y = x
  y[is.nan(y)] = 0
  return(y)
  }