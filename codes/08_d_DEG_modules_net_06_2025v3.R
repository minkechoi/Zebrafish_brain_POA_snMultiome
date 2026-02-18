
#module merge
#####
#load DEGs
load(file=paste0("./data/rda/DEGs_step7.rda"))

'

#loop for bPAC vs cont

bPAC_cont <- list()
for(i in names(list_contrasts2)){
  bPAC_cont[[i]] <-
    deseq_sc(
      m = dac_mat,
      d = sampletable,
      contrast_info = list_contrasts2[[i]],
      filter_by = "padj",
      p_threshold = 0.05,
      cell = i,
      plot_results = FALSE
    )
}


  list_contrasts,
  list_contrasts2,
  dac_mat,
  sampletable,
  postLD_preLD,
  bPAC_cont,
  '

#load tables for module info

library(readxl)

# read_excel reads both xls and xlsx files

top_genes_md=list()
for (i in libr) {
  load(file = paste0("data/",i,"/rda/danio_wgcna_all.rda"))
  
  gene_lists_by_column <- lapply(colnames(datKME), function(col) {
    genes <- rownames(datKME)[datKME[, col] > 0.7]
    val = datKME[genes,col]
    names(val)=genes
    val=sort(val, decreasing = TRUE)
    val <- head(names(val), 50)
    return(val)
  })
  names(gene_lists_by_column) <- colnames(datKME)
  
  top_genes_md[[i]]=gene_lists_by_column
}


merged_modules_sum=list()
for (i in libr) {
  module_sum=read_excel(paste0("./outputs/",i,"/ss_modules_table.xlsx"))
  module_sum$source=i
  merged_modules_sum[[i]]=module_sum
}

#data.frame
merged_modules_sum_tb=rbind(merged_modules_sum[[1]],merged_modules_sum[[2]],
                            merged_modules_sum[[3]],merged_modules_sum[[4]])


# Split and unnest the `celltypes` column
library(tidyr)
library(dplyr)

merged_modules_sum_tb <- merged_modules_sum_tb %>%
  separate_rows(celltypes, sep = ",")  # Split each celltypes entry into rows


# Group by celltypes and split
df_groups <- na.omit(merged_modules_sum_tb) %>%
  group_by(celltypes) %>%
  group_split()

# Extract celltype names in the correct order
celltype_names <- na.omit(merged_modules_sum_tb) %>%
  group_by(celltypes) %>%
  group_keys() %>%
  pull(celltypes)

#celltype_names=gsub("_","\\.",celltype_names)

# Assign names to the list
celltype_tables <- setNames(df_groups, celltype_names)

#select module genes by celltype
# Create an empty list to store output
extracted_genes <- list()

# Loop through each cell type group
for (celltype in names(celltype_tables)) {
  df <- celltype_tables[[celltype]]
  
  matched_genes=c()
  if (nrow(df)>0) {
    for (i in seq_len(nrow(df))) {
      mod <- df$newname[i]
      src <- df$source[i]
      
      # Extract corresponding genes from gene_table
      matching_genes <- top_genes_md[[src]][[mod]]
      
      # Store with a clear name like: galn__m01_cont_pre
      
      matched_genes=c(matched_genes,matching_genes)
    }
    key <- celltype
    extracted_genes[[key]] = unique(matched_genes)
  }
}

#for 45_sst1.1
extracted_genes[["45_sst1.1"]]=unique(c(extracted_genes[["45_sst1.1"]],"sst1.1"))

## Data Preparation
wgcna_counts=list()
##We prepare by loading the necessary data from our previous WGCNA and DEG:

#norm_count_matrix
for (group in libr) {
  #{r load_data}
  
  load(paste0("data/",group,"/rda/danio_counts.rda"))
  wgcna_counts[[group]]=danio_counts_norm
}
gc()

list_test_cluster=c("35.0_avp.crhb", "35.1_avp", "45_sst1.1")

DEG_and_module=c()
for (cluster in list_test_cluster) {
  ncluster=gsub("_","\\.",cluster)
  DEG_names= c(names(postLD_preLD)[grep(ncluster,names(postLD_preLD))],
               names(bPAC_cont)[grep(ncluster,names(bPAC_cont))])
  for (cl.key in DEG_names) {
    if (cl.key %in% names(postLD_preLD)) {
      ld_DEG=postLD_preLD[[cl.key]][["diffgenes"]]
      DEGs=c(DEGs,ld_DEG)
    }else{
      bPAC_DEG=bPAC_cont[[cl.key]][["diffgenes"]]
      DEGs=c(DEGs,bPAC_DEG)
    }
  }
  DEGs=na.omit(DEGs)
DEG_and_module=unique(c(DEGs,extracted_genes[[cluster]]))
}

mt0=as.data.frame(DEG_and_module)
rownames(mt0)=mt0[,1]
colnames(mt0)="gene"

mt1=as.data.frame(wgcna_counts[[1]])[DEG_and_module,]
mt1$gene=rownames(mt1)
mt1=mt1[,c("gene",list_test_cluster)]
colnames(mt1)=c("gene",paste0(list_test_cluster,"_",libr[1]))

mt2=as.data.frame(wgcna_counts[[2]])[DEG_and_module,]
mt2$gene=rownames(mt2)
mt2=mt2[,c("gene",list_test_cluster)]
colnames(mt2)=c("gene",paste0(list_test_cluster,"_",libr[2]))

mt3=as.data.frame(wgcna_counts[[3]])[DEG_and_module,]
mt3$gene=rownames(mt3)
mt3=mt3[,c("gene",list_test_cluster)]
colnames(mt3)=c("gene",paste0(list_test_cluster,"_",libr[3]))

mt4=as.data.frame(wgcna_counts[[4]])[DEG_and_module,]
mt4$gene=rownames(mt4)
mt4=mt4[,c("gene",list_test_cluster)]
colnames(mt4)=c("gene",paste0(list_test_cluster,"_",libr[4]))

library(dplyr)
library(purrr)

# Merge all by 'gene' using full_join
merged_df <- purrr::reduce(list(mt0,mt1,mt2,mt3,mt4), full_join, by = "gene")
merged_df=merged_df[-grep("NA", merged_df$gene),]
# Replace NAs with 0
merged_df[is.na(merged_df)] <- 0


# Set gene as rownames and remove column
rownames(merged_df) <- merged_df$gene
merged_df$gene <- NULL
merged_df=merged_df[rowSums(merged_df!=0) > 0, ]

write.csv(merged_df,paste0("./outputs/WGCNA_topQ_merged_mx.csv"),quote = F)

dir.create("./figures/cytoscape")

# Convert back to matrix (optional)
mt= as.matrix(t(scale(t(merged_df))))

a=Heatmap(
  name = "z-score_wgcna",
  mt,
  show_row_names=T,
  cluster_rows = T,km = 7,
  cluster_columns = T,
  col = colorRamp2(seq(0,1,length=6), rev(viridis_pastel)),
  #top_annotation = modules_ha
)
if (nrow(mt)>30) {
  pdf(paste0("./figures/cytoscape/wgcna_heatmap.pdf"),
       width = 5,height = 18)
  
  draw(a)
  
  dev.off()
}else{
  pdf(paste0("./figures/cytoscape/wgcna_heatmap.pdf"),
       width = 10,height = 10)
  
  draw(a)
  
  dev.off()
}


####

library(tidyverse)

# calculate correlation matrix
# the cor function works on columns so we have to transpose the matrix with the t() function
cor_matrix <- cor(t(merged_df))

# add gene ids and make a data frame
#colnames(cor_matrix) <- sig_data$Gene
cor_df <- tibble(
  Gene = colnames(cor_matrix),
  as.data.frame(cor_matrix)
)

# pivot the matrix to make it tidy
cor_long <- cor_df %>% 
  pivot_longer(cols = -Gene, names_to = 'Gene2', values_to = 'cor') %>% 
  # sort by Gene
  dplyr::arrange(Gene, Gene2) %>% 
  # remove where the gene is the same or higher
  # also only keep ones where cor > 0.8
  dplyr::filter(!(Gene == Gene2 | Gene > Gene2), cor > 0.8) %>% 
  # add a type column that is always "geneExprCor"
  dplyr::mutate(type = "geneExprCor") %>% 
  # reorder the columns
  dplyr::select(Gene, type, Gene2, cor)
# this is now the edges in the source -> target format

# write out network file
# change this filename if necessary
network_file <- paste0("wgcna_3545",'.sif')
write.table(dplyr::select(cor_long, Gene, type, Gene2), 
          paste0("./outputs/net_",network_file), col.names = F,row.names = F,quote = F,sep = "\t")

# get node info (DEG)
node_info <- merged_df
node_info$Gene=rownames(merged_df)

node_info_tb <- node_info %>%
  rowwise() %>%
  dplyr::mutate(
         sum = sum(c_across(-Gene), na.rm = TRUE),
         mean = mean(c_across(-Gene), na.rm = TRUE),
         median = median(c_across(-Gene), na.rm = TRUE),
         sd = sd(c_across(-Gene), na.rm = TRUE),
         max = max(c_across(-Gene), na.rm = TRUE)) %>%
  ungroup()
node_info_tb$sigDEG="no"
node_info_tb$sigDEG[which(node_info_tb$Gene %in% DEGs)]= "yes"
# write out nodes file
# change this filename if necessary
nodes_file <- paste0("wgcna_3545",'.txt')
write_tsv(node_info_tb, paste0("./outputs/node_",nodes_file))

# make edges df
unite(cor_long, 'edge', Gene, type, sep = ' (') %>% 
  unite(col = 'edge', edge, Gene2, sep = ') ') %>% 
  # write out edges file
  write_tsv(., paste0("./outputs/edges_",nodes_file)
           )

####

#heatmap for Fig3
cluster="45_sst1.1"
node_info_tb_45= read.csv("./outputs/45_sst1.1_WGCNA_topQ_merged_mx.csv",row.names = 1)
colnames(node_info_tb_45)=c("45_sst1.1_cont_pre","45_sst1.1_cont_post",
                              "45_sst1.1_bPAC_pre", "45_sst1.1_bPAC_post")
mt= as.matrix(t(scale(t(node_info_tb_45))))

a=Heatmap(
  name = cluster,
  mt,
  show_row_names=T,
  cluster_rows = T,
  cluster_columns = T,
  col = colorRamp2(seq(0,1.5,length=6), rev(viridis_pastel)),
  #top_annotation = modules_ha
)

pdf(paste0("./figures/WGCNA_",cluster,"_heatmap.pdf"), he = 15, wi = 5)
draw(a)
dev.off()


