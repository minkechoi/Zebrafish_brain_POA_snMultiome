chisq_and_posthoc <- function(tbl,method_ph="BH"){
  require(rstatix)
  require(chisq.posthoc.test)
  require(dplyr)
  require(tidyverse)
  require(tidyr)
  
  # chi squared
  chsq <- chisq_test(tbl)
  
  sample_cols <- ncol(tbl)
  
  # posthoc test
  ph <- chisq.posthoc.test(tbl, method = method_ph) %>% 
    pivot_longer(2+1:sample_cols,names_to = "category", values_to = "estimate") %>%
    dplyr::rename("Statistic"=Value, "Value"=estimate, "Cluster" = Dimension) %>%
    dplyr::select(Cluster,category, Statistic,Value)
  
  ## Aggregate at the cluster, category level(s) so pvalues and residuals are not in separate rows
  ph_tidy <- 
    left_join(
      ph[ph$Statistic=="p values",], 
      ph[ph$Statistic=="Residuals",],
      by = c("Cluster","category")
    ) %>%
    dplyr::rename("pval"=Value.x, "res"= Value.y) %>%
    dplyr::select(Cluster, category, pval, res) %>%
    arrange(Cluster, desc(category)) %>% 
    mutate(
      significance = cut(
        pval,
        breaks = c(0,0.001,0.01,0.05,1),
        labels = c("***","**","*"," "), include.lowest = T
      )
    )
  
  res <- 
    list(
     chi_squared_test = chsq,
     posthoc_test = ph,
     result = ph_tidy
     )
  
  return(res)
}