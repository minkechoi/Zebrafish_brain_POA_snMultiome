map_up <- function(x){
  a <- abs(x)
  a_max <- max(a)
  f <- magn_order(a_max)
  y <- x + 10**f
  return(y)
}