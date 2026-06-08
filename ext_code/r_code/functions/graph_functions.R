#' Comb graphs for nice plotting
top_edges_per_tf <- function(g,tf,top=3, mode = "out"){
  
  if(mode == "out"){
    tf_edges <- E(g)[.from(tf)]
  } else if(mode == "in"){
    tf_edges <- E(g)[.to(tf)]
  } else {
    Stop("incorrect mode. Check input parameters")
  }
  
  if (length(tf_edges) > top){
    tf_edges_top <- rev(order( E(g)[tf_edges]$weight ))[ 1:top ]    
  } else {
    tf_edges_top <- rev(order( E(g)[tf_edges]$weight ))
  }
  
  edges_top <- E(g)[which(E(g) %in% tf_edges)][tf_edges_top]
  
  y <- which(E(g) %in% edges_top)
  
  return(y)
}

#' Create a subgraph keeping the top interactions stemming from every TF
subgraph_with_top_edges_per_tf <- function(g, top = 3, delete_isolated = TRUE, mode = "out"){
  
  top_edges <- integer()
  
  if(mode %in% c("out","in")){
    for(tf in V(g)$name[igraph::degree(g,V(g), mode = "out") > 0]){ # > 0 == TFs
      e <- top_edges_per_tf(g=g, tf = tf, top = top, mode = mode)
      top_edges <- c(top_edges,e)
    }
  } else if(mode == "both"){
    # in
    e_in_top <- integer()
    for(tf in V(g)$name[igraph::degree(g,V(g), mode = "out") > 0]){ # > 0 == TFs
      e <- top_edges_per_tf(g=g, tf = tf, top = top, mode = "in")
      e_in_top <- c(e_in_top,e)
    }
    # out
    e_out_top <- integer()
    for(tf in V(g)$name[igraph::degree(g,V(g), mode = "out") > 0]){ # > 0 == TFs
      e <- top_edges_per_tf(g=g, tf = tf, top = top, mode = "out")
      e_out_top <- c(e_out_top,e)
    }
    
    top_edges <- c(e_in_top, e_out_top)
  }
  
  g_ <- subgraph.edges(graph = g, eids = E(g)[top_edges], delete.vertices = delete_isolated)
  
  return(g_)
}