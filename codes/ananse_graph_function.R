top_edges_per_tgtf <- function(g,tgtf,top=3){
  
  tg_edges <- E(g)[.to(tgtf)]
  tf_edges <- E(g)[.from(tgtf)]
  
  total_edge=c(tg_edges,tf_edges)
  
  if (length(tg_edges) > top){
    tg_edges_top <- rev(order( E(g)[tg_edges]$weight ))[ 1:top ]
  } else {
    tg_edges_top <- rev(order( E(g)[tg_edges]$weight ))
  }
  
  if (length(tf_edges) > top){
    tf_edges_top <- rev(order( E(g)[tf_edges]$weight ))[ 1:top ]
  } else {
    tf_edges_top <- rev(order( E(g)[tf_edges]$weight ))
  }
  
  total_edge_top=c(tg_edges_top,tf_edges_top)
  
  edges_top <- E(g)[which(E(g) %in% total_edge)][total_edge_top]
  
  y <- which(E(g) %in% edges_top)
  
  return(y)
  
}