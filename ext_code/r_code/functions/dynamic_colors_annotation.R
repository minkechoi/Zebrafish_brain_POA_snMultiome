#' @param df : a dataframe with columns as different categorisations, each with 
#' different colour palette and every unique event with an assigned color in that palette 
#' @param pals : a list of as many elements as columns in the data frame. Every element 
#' has the name of a given df column. Every element is a list of named colours. Colours 
#' are generated from a random categorical palette generator. Every name corresponds 
#' to a unique event in the given df column.
#' @param rand.seed : a numeric value to pass to the set.seed() command that will determine 
#' the random palettes generated.
dynamic_colors_annotation <- function(df, pals = NULL, rand.seed = 1234){
  library(colorspace)
  library(RColorBrewer)
  library(viridis)
  random_pals <- function(n){
    l_pals <- 
      list(
        function(x){return(brewer.pal(x,"Spectral"))},
        function(x){return(viridis(x))},
        function(x){return(cividis(x))},
        function(x){return(magma(x))},
        function(x){return(qualitative_hcl(x, "Warm"))},
        function(x){return(qualitative_hcl(x, "Cold"))},
        function(x){return(qualitative_hcl(x, "Set2"))},
        function(x){return(divergingx_hcl(x, "Earth"))}, 
        function(x){return(divergingx_hcl(x, "Fall"))}, 
        function(x){return(divergingx_hcl(x, "TealRose"))}, 
        function(x){return(divergingx_hcl(x, "Temps"))}, 
        function(x){return(divergingx_hcl(x, "RdYlBu"))}, 
        function(x){return(divergingx_hcl(x, "Roma"))}
      )
    return(l_pals[[sample(1:length(l_pals),1)]](n))
  }
  
  l <- list()
  
  if( is.null(pals) ){
    
    pals <-  list()
    
    set.seed(rand.seed)
    
    for(i in 1:ncol(df)){
      pals[[i]] <- random_pals(length(unique(df[,i])))
      names(pals)[i] <- colnames(df)[i]
    }
  }
  
  for(i in 1:ncol(df)){
    l[[i]] <- setNames(pals[[i]],unique(df[,i]))
    names(l)[i] <- colnames(df)[i]
  }
  
  return(l)
}
