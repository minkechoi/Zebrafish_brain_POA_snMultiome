top_edges_per_tf <- function(g,tf,top=3){
  
  tf_edges <- E(g)[.from(tf)]
  
  if (length(tf_edges) > top){
    tf_edges_top <- rev(order( E(g)[tf_edges]$weight ))[ 1:top ]    
  } else {
    tf_edges_top <- rev(order( E(g)[tf_edges]$weight ))
  }
  
  edges_top <- E(g)[which(E(g) %in% tf_edges)][tf_edges_top]
  
  y <- which(E(g) %in% edges_top)
  
  return(y)
  
}