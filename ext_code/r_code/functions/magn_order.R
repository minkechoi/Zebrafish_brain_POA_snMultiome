# from https://discuss.python.org/t/a-function-for-calculating-magnitude-order-of-a-number/18924
magn_order <- function(num){
  if (num == 0){
    return(0)
  }
  absnum <- abs(num)
  ord <- log10(absnum)
  res <- floor(ord)
  
  return(res)
}