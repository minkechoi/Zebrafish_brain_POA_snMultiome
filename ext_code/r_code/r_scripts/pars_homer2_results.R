homer_fd="./outputs/2025/homer_results/wgcna_homer2_atac_peak/known/"

#homer_fd="./outputs/2025/homer_results/wgcna_all/"
#homer_fd="./homer_result/fkbp5/"
known_HR_all=c(paste0(homer_fd,sort(list.files(homer_fd))))
#known_HR_all=c(paste0(homer_fd,sort(list.files(homer_fd)),
#         "/knownResults.txt"))
modules=sort(list.files(homer_fd))
#modules= str_split(string = sort(list.files(homer_fd)),
#                   pattern = "_") %>% sapply(function(x) x[1])
#modules= str_split(string = sort(list.files(homer_fd)),
#                   pattern = "\\.") %>% sapply(function(x) x[1])
#load know motif files
i=1
tsv_motifs_all=vroom::vroom(known_HR_all[i])
tsv_motifs_all$target_peakset_name = modules[i]
tsv_motifs_all=tsv_motifs_all %>% relocate(target_peakset_name, .before = everything())
colnames(tsv_motifs_all) <- 
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
for (i in 2:length(known_HR_all)) {
  tsv_motifs=vroom::vroom(known_HR_all[i])
  tsv_motifs$target_peakset_name = modules[i]
  tsv_motifs=tsv_motifs %>% relocate(target_peakset_name, .before = everything())
  
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
  
  tsv_motifs_all=rbind(tsv_motifs_all,tsv_motifs)
}

categ_regex1=""
categ_regex2=""
qval_thresh = 0.1 
max_logqval = 10

homer_trim=parse_homer2_output_table(tsv_motifs_all)
