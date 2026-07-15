# =============================================================================
# ananse_graph_function.R
# -----------------------------------------------------------------------------
# Purpose : Helper for ANANSE gene-regulatory-network graphs (igraph). Given a
#           node `tgtf` (a TF or target gene), return the indices of its highest-
#           weight edges: the top incoming edges (regulators of the node) and the
#           top outgoing edges (targets it regulates).
# Args    : g     - an igraph GRN graph with edge attribute `weight`
#           tgtf  - node of interest
#           top   - max number of edges to keep in each direction (default 3)
# Returns : integer indices into E(g) of the selected top edges.
# =============================================================================

top_edges_per_tgtf <- function(g,tgtf,top=3){

  tg_edges <- E(g)[.to(tgtf)]     # edges pointing INTO tgtf (its regulators)
  tf_edges <- E(g)[.from(tgtf)]   # edges going OUT of tgtf (its targets)

  total_edge=c(tg_edges,tf_edges) # all edges incident to the node

  # Keep the `top` strongest incoming edges by weight (or all if fewer exist)
  if (length(tg_edges) > top){
    tg_edges_top <- rev(order( E(g)[tg_edges]$weight ))[ 1:top ]
  } else {
    tg_edges_top <- rev(order( E(g)[tg_edges]$weight ))
  }

  # Keep the `top` strongest outgoing edges by weight (or all if fewer exist)
  if (length(tf_edges) > top){
    tf_edges_top <- rev(order( E(g)[tf_edges]$weight ))[ 1:top ]
  } else {
    tf_edges_top <- rev(order( E(g)[tf_edges]$weight ))
  }

  total_edge_top=c(tg_edges_top,tf_edges_top) # combined selection (positions within total_edge)

  edges_top <- E(g)[which(E(g) %in% total_edge)][total_edge_top] # map back to actual edges

  y <- which(E(g) %in% edges_top)  # final edge indices in the full graph

  return(y)

}
