

parse_homer_output_table <- function(tsv_motifs,categ_regex1 = "", categ_regex2 = "", qval_thresh = 0.1, max_logqval = 10){
  library(tidyverse)
  library(stringr)
  
  # Read the file and filter out rows with "Consensus" in the "Consensus" column
  tsv_motifs <- tsv_motifs %>%
    dplyr::filter(!str_detect(Consensus, "Consensus"))
  
  # Rename columns and clean up the names of the clusters
  colnames(tsv_motifs) <- 
    c(
      "target_peakset_name",
      "motif",
      "Consensus",
      "pvalue",
      "logpval",
      "qval",
      "no_target_seqs_with_motif",
      "pct_target_seqs_with_motif",
      "no_bg_seqs_with_motif",
      "pct_bg_seqs_with_motif"
    )
  
  ## Filter rows based on a threshold value and calculate new columns from existing columns
  tsv_motifs$pct_target_seqs_with_motif <- gsub("%","",tsv_motifs$pct_target_seqs_with_motif) #i had these using a smartsy of mutate+across+str_replace_all but it was not working so I reverted back to base.
  tsv_motifs$pct_bg_seqs_with_motif <- gsub("%","",tsv_motifs$pct_bg_seqs_with_motif)
  
  tsv_motifs <- tsv_motifs %>%
    mutate_at(vars(pvalue:pct_bg_seqs_with_motif), as.numeric) %>%
    dplyr::filter(qval < qval_thresh ) %>%
    mutate(
      logpval = -1 * logpval,
      logqval = ifelse(qval == 0, max_logqval, -log(qval)),
      motif = gsub("/.*","",motif),
      class = gsub("(.*\\((.*)\\).*)", "\\2",motif) %>% str_remove_all("[()]"),
      target_peakset_name = gsub(categ_regex1,categ_regex2,target_peakset_name)
    )
  
  
  total_num_peaks <- function(x){
    num <- tsv_motifs$no_bg_seqs_with_motif[tsv_motifs$target_peakset_name == x ]
    pct <- tsv_motifs$pct_target_seqs_with_motif[tsv_motifs$target_peakset_name == x ]
    abc <- num / pct
    rs <- round(mean(abc[!is.infinite(abc) & !(is.nan(abc))])*100)
    return(rs)
  }
  
  tgt_num_peaks <- function(x){
    num <- tsv_motifs$no_target_seqs_with_motif[tsv_motifs$target_peakset_name == x ]
    pct <- tsv_motifs$pct_target_seqs_with_motif[tsv_motifs$target_peakset_name == x ]
    abc <- num / pct
    rs <- round(mean(abc[!is.infinite(abc) & !(is.nan(abc))])*100)
    return(rs)
  }
  
  num_background_peaks <- sapply(unique(tsv_motifs$target_peakset_name),total_num_peaks)
  num_tgt_peaks <- sapply(unique(tsv_motifs$target_peakset_name),tgt_num_peaks)
  
  tsv_motifs$num_total_tgt_seqs <- sapply(tsv_motifs$target_peakset_name,function(x){return(num_tgt_peaks[names(num_tgt_peaks) == x])})
  
  tsv_motifs$num_total_bg_seqs <- sapply(tsv_motifs$target_peakset_name,function(x){return(num_background_peaks[names(num_background_peaks) == x])})
  
  tsv_motifs$size_pct <- as.numeric(cut(tsv_motifs$pct_target_seqs_with_motif,breaks = c(0,10,20,40,Inf), labels = c(1,2,3,4)))
  
  # USE A SMART DICTIONARY FUNCTION FOR THIS
  #merging together a bunch of motifs to prevent over-representation
  tsv_motifs$motif_minim <- tidyup_motifnames(tsv_motifs$motif)
  
  return(tsv_motifs)
}



parse_homer2_output_table <- function(tsv_motifs,categ_regex1 = "", categ_regex2 = "", qval_thresh = 0.1, max_logqval = 10){
  library(tidyverse)
  library(stringr)
  
  # Read the file and filter out rows with "Consensus" in the "Consensus" column
  tsv_motifs <- tsv_motifs %>%
    dplyr::filter(!str_detect(Consensus, "Consensus"))
  
  # Rename columns and clean up the names of the clusters
  colnames(tsv_motifs) <- 
    c(
      "target_peakset_name",
      "motif",
      "Consensus",
      "pvalue",
      "logpval",
      "log2Enrich",
      "no_target_seqs_with_motif",
      "pct_target_seqs_with_motif",
      "no_bg_seqs_with_motif",
      "pct_bg_seqs_with_motif"
    )
  
  ## Filter rows based on a threshold value and calculate new columns from existing columns
  tsv_motifs$pct_target_seqs_with_motif <- gsub("%","",tsv_motifs$pct_target_seqs_with_motif) #i had these using a smartsy of mutate+across+str_replace_all but it was not working so I reverted back to base.
  tsv_motifs$pct_bg_seqs_with_motif <- gsub("%","",tsv_motifs$pct_bg_seqs_with_motif)
  
  tsv_motifs <- tsv_motifs %>%
    mutate_at(vars(pvalue:pct_bg_seqs_with_motif), as.numeric) %>%
    dplyr::filter(2^-log2Enrich < qval_thresh ) %>%
    mutate(
      logpval = -1 * logpval,
      logqval = ifelse(2^-log2Enrich == 0, max_logqval, -log(2^-log2Enrich)),
      motif = gsub("/.*","",motif),
      class = gsub("(.*\\((.*)\\).*)", "\\2",motif) %>% str_remove_all("[()]"),
      target_peakset_name = gsub(categ_regex1,categ_regex2,target_peakset_name)
    )
  
  
  total_num_peaks <- function(x){
    num <- tsv_motifs$no_bg_seqs_with_motif[tsv_motifs$target_peakset_name == x ]
    pct <- tsv_motifs$pct_target_seqs_with_motif[tsv_motifs$target_peakset_name == x ]
    abc <- num / pct
    rs <- round(mean(abc[!is.infinite(abc) & !(is.nan(abc))])*100)
    return(rs)
  }
  
  tgt_num_peaks <- function(x){
    num <- tsv_motifs$no_target_seqs_with_motif[tsv_motifs$target_peakset_name == x ]
    pct <- tsv_motifs$pct_target_seqs_with_motif[tsv_motifs$target_peakset_name == x ]
    abc <- num / pct
    rs <- round(mean(abc[!is.infinite(abc) & !(is.nan(abc))])*100)
    return(rs)
  }
  
  num_background_peaks <- sapply(unique(tsv_motifs$target_peakset_name),total_num_peaks)
  num_tgt_peaks <- sapply(unique(tsv_motifs$target_peakset_name),tgt_num_peaks)
  
  tsv_motifs$num_total_tgt_seqs <- sapply(tsv_motifs$target_peakset_name,function(x){return(num_tgt_peaks[names(num_tgt_peaks) == x])})
  
  tsv_motifs$num_total_bg_seqs <- sapply(tsv_motifs$target_peakset_name,function(x){return(num_background_peaks[names(num_background_peaks) == x])})
  
  tsv_motifs$size_pct <- as.numeric(cut(tsv_motifs$pct_target_seqs_with_motif,breaks = c(0,10,20,40,Inf), labels = c(1,2,3,4)))
  
  # USE A SMART DICTIONARY FUNCTION FOR THIS
  #merging together a bunch of motifs to prevent over-representation
  tsv_motifs$motif_minim <- tidyup_motifnames(tsv_motifs$motif)
  
  return(tsv_motifs)
}


tidyup_motifnames <- function(x){
  a <- str_split(string =x,
                 pattern = "/") %>% sapply(function(b) b[1])
  
  return(a)
}

aggregate_motifs_table <- function(tsv, fun = max){
  # took the idea from this R API for Homer : https://robertamezquita.github.io/marge/articles/marge-workflow.html
  aggreg <- 
    aggregate(
      tsv[,c(8,11,15)],
      by = list(
        class = tsv$class,
        category = tsv[[1]]),
      FUN = fun
    ) 
  
  colnames(aggreg) <- 
    c("class","category","pct","logqval","size_pct")
  
  return(aggreg)
}