smed_cisreg_all_tss <- read.table("/mnt/sda/alberto/genomes/Smed/Rink/schMedS3/schMedS3_h1_ENCODE_hybrid_agat_hconf_TSS.bed6.bed", header = FALSE)
smed_rosetta <- read.delim2("/mnt/sda/alberto/projects/smed_rink_gene_annot/20230804_Smed_Rink_Simplified_Annotation_Table.tsv", header = TRUE)

smed_cisreg_longestIso_tss <- 
  smed_cisreg_all_tss[smed_cisreg_all_tss$V4 %in% smed_rosetta$longest_isoform,]

smed_cisreg_tss <- smed_cisreg_longestIso_tss
smed_cisreg_tss$V4 <- 
  translate_ids(x = smed_cisreg_tss$V4, dict = smed_rosetta[,c(2,1)])

write.table(
  smed_cisreg_tss,
  file = "/mnt/sda/alberto/genomes/Smed/Rink/schMedS3/schMedS3_h1_ENCODE_hybrid_agat_hconf_longestIso_usingGeneID_TSS.bed6.bed",
  sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE
  )
