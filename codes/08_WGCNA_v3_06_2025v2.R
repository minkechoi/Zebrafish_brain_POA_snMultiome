
# =============================================================================
# 08_WGCNA_v3_06_2025v2.R
# -----------------------------------------------------------------------------
# Purpose : Weighted gene co-expression network analysis (WGCNA) driver. For each
#           library it runs the pseudobulk-prep + WGCNA pipeline (functions in
#           08_b), estimates the soft-threshold power, detects co-expression
#           modules, relates modules to cell types and to the stress signatures
#           (AUCell/UCell), tests module TF connectivity, and runs GO + Signac
#           motif enrichment. Also defines run_merged_WGCNA_GO() for merged runs.
# Depends : 08_a (pseudobulk prep), 08_b (WGCNA functions), 08_c (graph functions).
# Inputs  : per-library danio_counts.rda / danio_wgcna_all.rda (from 08_a/08_b).
# Outputs : module tables, module-trait heatmaps, GO/motif results, and figures
#           under ./outputs/<lib>/ and ./figures/<lib>/WGCNA/.
# Sections: setup -> WGCNA (power estimation -> modules) -> visualization ->
#           TF connectivity / WGCNA graph -> UCell module scoring -> GO -> motif.
# Author  : Min K Choi, m.choi@exeter.ac.uk
# =============================================================================
#Author: Min K Choi, m.choi@exter.ac.uk
#DEG analysis

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc()

##WGCNA analysis

## Load libraries

library(plyr)
library(dplyr)
library(ggplot2)
library(Seurat)
library(colorspace)
library(SCP)
#set the data and annotations 
set.seed(1234)
options(future.globals.maxSize = 120000 * 1024^2)

r.variable=4000
vs="06_2025_v4"
umap="wnn.umap"

#set working dir
if(getwd()!=paste0("D:/projects/scMultiome_oxt/",vs)){
  setwd(paste0("./",vs))
}

#functions adopted and modified from Alberto Perez-Posada @apposada
source("./ext_code/r_code/functions/sourcefolder.R")
sourceFolder("./ext_code/r_code/functions/")

#other temp function
'%!in%' <- function(x,y)!('%in%'(x,y))

##WGCNA analysis
## Load Data
ss1=readRDS(paste0("./data/rds/Step7_var",r.variable,".rds"))
DefaultAssay(ss1)="RNA"
libr=unique(ss1$orig.ident)
type_table_m.v2= read.csv(paste0("./outputs/cell_type/cell_type_table_m_modified.csv"), header = T, row.names = 1)
type_table_m.v3= type_table_m.v2[- c(grep(pattern = "UnD", type_table_m.v2$numbered)),]

libr=unique(ss1$orig.ident)
ctype_path=paste0("./outputs/cell_type/ss4000_ctypes.csv")
matrix_type= "danio_counts_norm" # choose use cw or not: danio_counts_norm_cw or danio_counts_norm
hclust.method="ward.D2" # choose method: "ward.D","ward.D2", "complete", "average", ...

source("./codes/08_b_WGCNA_functions.R")

for (i in libr[1:4]) {
  print(i)
  wgcna_prep(ss1,i,ctype_path,matrix_type,hclust.method)
}


# WGCNA  ------------------------------------------------------------------


## About

#This markdown showcases the use of WGCNA (Langfelder & Horvath, 2008) to identify modules of genes that are regulated in a similar manner across clusters of cell types, using pseudo-bulk data from gene counts aggregated at the Leiden cluster level.
#This markdown contains the main code chunks necessary for data transformation and running of WGCNA. For more information and a more detailed documentation, please refer to this script from the original documentation:
#  https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-NetworkConstruction-man.R

## Loading Necessary Packages

#{r load_packages, warning = FALSE, message=FALSE}
library(data.table)
library(reshape2)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(viridis)
library(colorspace)
library(WGCNA)
library(topGO)
library(plyr)
library(xlsx)

## Loading necessary code (CHANGE PATHS)

## Data Preparation

#We prepare by loading the necessary data from our previous markdowns:


#We load our table of schmidtea counts as input for the wgcna analysis. 
#We will select variable genes to enter the analysis of WGCNA, with a CV > 1. Gene expression data will also be scaled by centering around the mean (z-score). datExpr is the data frame in the format that WGCNA likes.

#optimal cv calulation
for (i in libr) {
load(paste0("./data/",i,"/rda/danio_counts.rda"))
 
##{r, message=FALSE}
danio_cv <- apply(danio_counts_norm_cw,1,function(x){sd(x)/mean(x)})
a=plot(density(danio_cv), main = "CV")
a+abline(v=quantile(danio_cv), col = brewer.pal(6,"Spectral")[-4], lty = 2, lwd = rep(2))+
  a+abline(v=c(1,1.25,1.5,2), col = c("#30a958","#a0d16d","#e99f4e","#d34646"), lwd = 1.5)

pdf(paste0("./figures/",i,"/WGCNA/danio_cv.pdf"),
     width = 15,height = 8)
plot(density(danio_cv), main = "CV")+abline(v=quantile(danio_cv), col = brewer.pal(6,"Spectral")[-4], lty = 2, lwd = rep(2))+abline(v=c(1,1.25,1.5,2), col = c("#30a958","#a0d16d","#e99f4e","#d34646"), lwd = 1.5)

dev.off()
}

i=1
#softpower calulation
for (i in 1:4) {
load(paste0("./data/",libr[i],"/rda/danio_counts.rda"))
danio_cv <- apply(danio_counts_norm_cw,1,function(x){sd(x)/mean(x)})
print(quantile(danio_cv)[3])
thresh_cv <- 1 #quantile(danio_cv)[3], 1.25
filt_cv <- which(danio_cv >= thresh_cv)
input_wgcna <- danio_counts_norm_cw[filt_cv,]

# Defining the datExpr object
datExpr = as.data.frame(scale(t(input_wgcna)))
rownames(datExpr) <- colnames(danio_counts_norm_cw)


### Power estimation : pick the soft-threshold power giving a scale-free topology

#We generate a set of numbers to run the analysis of scale free topology.

#{r}
# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
powers


#We proceed to run the analysis of scale free topology for multiple soft thresholding powers. The aim is to help the user pick an appropriate soft-thresholding power for network construction.

#{r power_estimation, message = FALSE}
# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


#Browsing the fitting of the topology model we can pick our Beta soft thresholding power of choice, as it provides a high R^2 fit and it also lays at the turning point between a lot of connectivity and very low connectivity.

#{r , fig.width=8, fig.height=5}
plot_scalefreetopology_pretty(sft)

pdf(paste0("./figures/",libr[i],"/WGCNA/sft_thresholding.pdf"),wi = 8, he = 5)
plot_scalefreetopology_pretty(sft)
dev.off()

save(
  danio_ctypes,
  danio_counts,
  danio_psbulk_cellweights,
  danio_counts_norm,
  danio_counts_norm_cw,
  danio_genecolor,
  input_wgcna,
  datExpr,
  file = paste0("./data/",libr[i],"/rda/danio_counts.rda")
)

}

#With this we decide the number to which set the soft Power parameter:

#WGCNA run
hclust.method="ward.D2" # choose method: "ward.D","ward.D2", "complete", "average", ...

sps=c(4,6,6,9)
for (i in 1:4) {
run_WGCNA(ss1,sp=sps[i],hclust.md=hclust.method,min.md.size=30,group=libr[i])
}


#WGCNA TF connectivity
mkss= 0.65
for (i in 1:4) {
  run_WGCNA_TF(obj=ss1,mks=mkss,hclust.md=hclust.method,min.md.size=30,group=libr[i])
}

#WGCNA module list export for homer
for (i in 1:4) {
  module_out(group=libr[i])
}


#WGCNA module GO enrichment test
for (i in 1:4) {
run_WGCNA_GO(group=libr[i])
}




####
#module_compare
#libr

library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(ggplot2)
library(effects)

library(readxl)
cont_preLD_md=read_excel(paste0("./outputs/",libr[1],"/danio_wgcna_id_module.xlsx"))
cont_postLD_md=read_excel(paste0("./outputs/",libr[2],"/danio_wgcna_id_module.xlsx"))
bPAC_preLD_md=read_excel(paste0("./outputs/",libr[3],"/danio_wgcna_id_module.xlsx"))
bPAC_postLD_md=read_excel(paste0("./outputs/",libr[4],"/danio_wgcna_id_module.xlsx"))

cont_preLD_md$source_md=paste0(libr[1],"_",cont_preLD_md$module)
cont_postLD_md$source_md=paste0(libr[2],"_",cont_postLD_md$module)
bPAC_preLD_md$source_md=paste0(libr[3],"_",bPAC_preLD_md$module)
bPAC_postLD_md$source_md=paste0(libr[4],"_",bPAC_postLD_md$module)

preLD_md=rbind(cont_preLD_md,bPAC_preLD_md)[,-2]
postLD_md=rbind(cont_postLD_md,bPAC_postLD_md)[,-2]

colnames(preLD_md)=c("id","module")
colnames(postLD_md)=c("id","module")

comp_LD=comparemodules(preLD_md,postLD_md)


###visualization


col_danio_tfs_expr <- colorRamp2(
  c(0:3),
  colorRampPalette(c("#f1f5ff","#b1b8c4","#2c2b46"))(4)
)

#module info

cont_preLD_md_link=read_excel(paste0("./outputs/",libr[1],"/ss_modules_table.xlsx"))
cont_postLD_md_link=read_excel(paste0("./outputs/",libr[2],"/ss_modules_table.xlsx"))
bPAC_preLD_md_link=read_excel(paste0("./outputs/",libr[3],"/ss_modules_table.xlsx"))
bPAC_postLD_md_link=read_excel(paste0("./outputs/",libr[4],"/ss_modules_table.xlsx"))

cont_preLD_md_link$source_md=paste0(libr[1],"_",cont_preLD_md_link$newname )
cont_postLD_md_link$source_md=paste0(libr[2],"_",cont_postLD_md_link$newname )
bPAC_preLD_md_link$source_md=paste0(libr[3],"_",bPAC_preLD_md_link$newname )
bPAC_postLD_md_link$source_md=paste0(libr[4],"_",bPAC_postLD_md_link$newname )

md_link=rbind(cont_preLD_md_link,cont_postLD_md_link,
              bPAC_preLD_md_link,bPAC_postLD_md_link)

hyge_mtxs=-log10(comp_LD$hypgeom)
danio_ctypes <- read.csv(ctype_path,row.names = 1)
md_link$single_ctype=translate_ids(md_link$cell_color,dict = danio_ctypes[,c(3,1)])
md_link$comp_name=paste(md_link$source_md,md_link$celltypes,sep = ":")


ctypes_rowAnno <-
  rowAnnotation(
    cluster = rownames(hyge_mtxs),
    col = list( cluster = setNames(md_link$cell_color ,md_link$source_md)),
    show_legend = F,show_annotation_name = F
  )

clu_ha = HeatmapAnnotation(
  name = "m_types",
  cluster = factor(colnames(hyge_mtxs)),
  col = list(cluster = setNames(md_link$cell_color ,md_link$source_md)),
  show_legend = F, 
  show_annotation_name = F)


  mt=hyge_mtxs
  h <- Heatmap(
    name=paste0("z-score"),
    t(scale(t(mt))),
    col = col_danio_tfs_expr,
    show_row_names = T,
    show_column_names = T,
    cluster_rows=T,clustering_method_rows = "complete", #"complete", "ward.D2","average
    cluster_columns=T, clustering_method_columns = "complete",
    top_annotation = clu_ha,
    bottom_annotation = clu_ha,
    right_annotation = ctypes_rowAnno,
    left_annotation = ctypes_rowAnno,
    #row_km =4,
    #column_km =5,
    cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height, 
              gp = gpar(col = "white", fill = NA, lty = 1))
    },
    row_title="post-LD vs. pre-LD"
  )
  draw(h)
  
  lgd=Legend(at=unique(md_link$single_ctype),title = "ctypes", 
             legend_gp = gpar(fill = setNames(unique(md_link$cell_color) ,
                                              unique(md_link$single_ctype))),
             #nrow=17,
             ncol = 1
             )


  pdf(paste0("./figures/WGCNA_module_compare_hypergeo1.pdf"),
       width = 50,height = 40)
  draw(h)

  dev.off()
  
  pdf(paste0("./figures/WGCNA_module_compare_hypergeo_cell_legend.pdf"),
       width = 3,height = 15)
  
  draw(lgd)
  
  dev.off()
  
#####
  
  hyge_mtxs=-log10(comp_LD$hypgeom)
  danio_ctypes <- read.csv(ctype_path,row.names = 1)
  md_link$single_ctype=translate_ids(md_link$cell_color,dict = danio_ctypes[,c(3,1)])
  md_link$comp_name=paste(md_link$source_md,md_link$celltypes,sep = ":")
  
  rn_mt=md_link[which(md_link$source_md %in% rownames(hyge_mtxs)),c(9,10,11)]
  cn_mt=md_link[which(md_link$source_md %in% colnames(hyge_mtxs)),c(9,10,11)]
  rownames(rn_mt)=rn_mt$source_md
  rownames(cn_mt)=cn_mt$source_md
  
  hyge_mtxs2=hyge_mtxs[rownames(rn_mt),rownames(cn_mt)]
  rownames(hyge_mtxs2)=rn_mt$comp_name
  colnames(hyge_mtxs2)=cn_mt$comp_name
  
  ctypes_rowAnno <-
    rowAnnotation(
      cluster = rownames(hyge_mtxs2),
      col = list( cluster = setNames(md_link$cell_color ,md_link$comp_name)),
      show_legend = F,show_annotation_name = F
    )
  
  clu_ha = HeatmapAnnotation(
    name = "m_types",
    cluster = factor(colnames(hyge_mtxs2)),
    col = list(cluster = setNames(md_link$cell_color ,md_link$comp_name)),
    show_legend = F, 
    show_annotation_name = F)
  
  
  mt=hyge_mtxs2
  h <- Heatmap(
    name=paste0("z-score"),
    t(scale(t(mt))),
    col = col_danio_tfs_expr,
    show_row_names = T,
    show_column_names = T,
    cluster_rows=T,clustering_method_rows = "complete", #"complete", "ward.D2","average
    cluster_columns=T, clustering_method_columns = "complete",
    top_annotation = clu_ha,
    bottom_annotation = clu_ha,
    right_annotation = ctypes_rowAnno,
    left_annotation = ctypes_rowAnno,
    #row_km =4,
    #column_km =5,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.rect(x = x, y = y, width = width, height = height, 
                gp = gpar(col = "white", fill = NA, lty = 1))
    },
    row_title="post-LD vs. pre-LD"
  )
  draw(h)
  
  lgd=Legend(at=unique(md_link$single_ctype),title = "ctypes", 
             legend_gp = gpar(fill = setNames(unique(md_link$cell_color) ,
                                              unique(md_link$single_ctype))),
             #nrow=17,
             ncol = 1
  )
  
  
  pdf(paste0("./figures/WGCNA_module_compare_hypergeo.pdf"),
       width = 50,height = 40)
  par(mar=c(100,2,2,100))
  draw(h)
  
  dev.off()
#####
#motif

#load homer results, Known motifs
motif_results=list()
for (i in 1:4) {
  paths=paste0("./outputs/",libr[i],"/results/")

  mtf=pars_homer(pth=paths)
  motif_results[[libr[i]]]=mtf
}


head(motif_results)



#set threshold
categ_regex1=""
categ_regex2=""
qval_thresh = 0.1 
max_logqval = 10


trimed_motif_results=list()
for (i in 1:4) {
  trimed_motif_results[[i]]=parse_homer_output_table(motif_results[[i]])
}

head(trimed_motif_results)

#motif_dotplots
for (i in 1:4) {
  motif_enrichdot(trimed_motif_results[[i]],libr[[i]])
}

####TF connectivity : rank transcription factors by intramodular connectivity (hub TFs)

for (i in libr) {
TF_primed(i)
}

####WGCNA graph
source("./codes/08_c_WGCNA_graph_functions.R")

for (i in libr) {
  wgcna_graph(i)
}

#module merge
#####
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
merged_modules_sum_tb=as.data.frame(merged_modules_sum_tb)

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

save(extracted_genes,
     celltype_tables,
     file = paste0("./data/rda/wgcna_genes.rda"))

library(openxlsx)
write.xlsx(extracted_genes, file = "./outputs/merged_wgcna_genes.xlsx")

celltype_tables_all=do.call(rbind, celltype_tables)
write.csv(celltype_tables_all,"./outputs/merged_wgcna_celltype_tables.csv",row.names = F)


####UCell : relate WGCNA modules to the stress signatures via UCell module scores
library(UCell)

#score calculation
#Ucell
DefaultAssay(ss1)="RNA"
for (i in names(extracted_genes)) {

signatures= extracted_genes[i]
ss1 <- AddModuleScore_UCell(ss1, 
                            features=signatures, name="_module_UCS")

#ss1 <- SmoothKNN(ss1,
#                 signature.names = paste0(i,"_module_UCS"),
#                 reduction="harmony",k=20, suffix = "_sm")

}

signature_plot=list()
for (i in names(extracted_genes)) {
  
signature_plot[[i]]= FeatureDimPlot(
    srt = ss1, features = paste0(i,"_module_UCS"),
    assay = "RNA",
    #label_repel = T,label_repulsion = 50,
    pt.size = 0.7,palcolor =c("lightgray","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
    seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, 
    #add_density = T,
    reduction = umap, theme_use = "theme_blank"
  )
}

signature_plot.split=list()
for (i in names(extracted_genes)) {
  
  signature_plot.split[[i]]= FeatureDimPlot(
    srt = ss1, features = paste0(i,"_module_UCS"),split.by = "orig.ident",
    assay = "RNA",nrow = 1,
    #label_repel = T,label_repulsion = 50,
    #pt.size = 0.7,
    palcolor =c("lightgray","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026"),
    seed = 0, compare_features = T, label =F, label_repel = T,label_insitu = TRUE, 
    #add_density = T,
    reduction = umap, theme_use = "theme_blank"
  )
}

pdf(paste0("./figures/module_signatureplot.pdf"),
     width = 25,height = 40)
plot_grid(plotlist = signature_plot,ncol = 5)
dev.off()


pdf(paste0("./figures/module_signatureplot.split.pdf"),
     width = 20,height = 45)
plot_grid(plotlist = signature_plot.split,ncol = 2)
dev.off()

#signiture heatmap
libs=unique(ss1$orig.ident)
load(file=paste0("./data/rda/hcluster.rda"))

auc_matrix_tb= ss1[[names(ss1@meta.data)[grepl("module_UCS",names(ss1@meta.data))]]]
texa_order_num=as.character(sapply(str_split(taxa_order[- grep("UnD",taxa_order)],"_"), function (x){x[[1]]}))
texa_order_num=gsub("\\.0","",texa_order_num)

#texa_order_num=factor(texa_order_num,levels = texa_order_num)
#auc
# Get column names to summarize (exclude the first two: merged_sub_numb and orig.ident)
value_cols <- setdiff(names(auc_matrix_tb), c("merged_sub_numb", "orig.ident"))

# Summarise using summarise_at (compatible with older versions of dplyr)
auc_matrix_tb$merged_sub_numb=factor(ss1$merged_sub_numb,levels =texa_order_num )
auc_matrix_tb$orig.ident=factor(ss1$orig.ident)

auc_sum_matrix_tb <- auc_matrix_tb %>%arrange(merged_sub_numb)%>%
  group_by(merged_sub_numb, orig.ident) %>%
  summarise_at(vars(value_cols), ~mean(., na.rm = TRUE)) %>%
  ungroup()


auc_sum_matrix=t(auc_sum_matrix_tb[,-c(1:2)])
colnames(auc_sum_matrix)=paste0(auc_sum_matrix_tb$merged_sub_numb,"_",auc_sum_matrix_tb$orig.ident)
rownames(auc_sum_matrix)=sapply(str_split(rownames(auc_sum_matrix),"_"),FUN = function(x){x[[1]]})
rownames(auc_sum_matrix)=as.numeric(rownames(auc_sum_matrix))
texa_reorder_num=texa_order_num[which(texa_order_num %in% rownames(auc_sum_matrix))]
auc_sum_matrix=auc_sum_matrix[c(texa_reorder_num),]

#visualization
#auc_heatmaps
#####
type_table_m.v2= read.csv(paste0("./outputs/cell_type/cell_type_table_m_modified.csv"), header = T, row.names = 1)
type_table_m.v3= type_table_m.v2[- c(grep(pattern = "UnD", type_table_m.v2$numbered)),]
type_table_m.v3$cl_numb=sapply(str_split(type_table_m.v3$numbered,"_"),function(x){x[[1]]})
type_table_m.v3$cl_numb=gsub("\\.0","",type_table_m.v3$cl_numb)

mt=t(scale(t(as.matrix(auc_sum_matrix))))

clu_ha = HeatmapAnnotation(
  name = "library",
  cluster = factor(colnames(mt), levels = unique(colnames(mt))),
  col = list(cluster = setNames(rep(rev(c("#d53e4f","#fee08b","#abdda4","#3288bd")),length(unique(colnames(mt)))/4),
                                unique(colnames(mt)))),
  show_legend = F, show_annotation_name = F
)

clu_ha2 = HeatmapAnnotation(
  name = "library",
  cluster = factor(colnames(mt), levels = unique(colnames(mt))),
  col = list(cluster = setNames(rep(type_table_m.v3$clcol,each=4),
                                as.vector((sapply(type_table_m.v3$cl_numb, function(x) paste0(x,"_",libr))))
                                )),
  show_legend = F, show_annotation_name = F
)

type.table=type_table_m.v3[which(as.character(type_table_m.v3$cl_numb) %in% rownames(mt)),]
  
ctypes_rowAnno <-
  rowAnnotation(
    cluster = rownames(mt),
    col = list( cluster = setNames(type.table$clcol,type.table$cl_numb) ),
    show_legend = F, show_annotation_name = F
  )

auc_signature <- ComplexHeatmap::Heatmap(
  name = "auc_signature",
  mt, #danio_wg_module_viz[,-modulecolumn], 
  cluster_rows= F,
  show_row_names = T,
  row_names_side = "left",
  show_row_dend = F,
  cluster_columns = F,
  column_title = "module top signature score",
  column_title_side = "bottom",
  show_column_names = F,
  column_names_side = "bottom",
  column_names_rot = 90,
  column_split =  factor(auc_sum_matrix_tb$merged_sub_numb),
  row_title_side = "left",
  row_title_rot = 0,
  border = T,border_gp = gpar(col = "grey", lty = 1),
  top_annotation = clu_ha2,
  left_annotation = ctypes_rowAnno,
  bottom_annotation = clu_ha,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x = x, y = y, width = width, height = height, 
              gp = gpar(col = "grey", fill = NA, lty = 1))
  },
  use_raster = FALSE,
  col = colorRampPalette(c("white","#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026"))(6)
)
draw(auc_signature)

pdf(paste0("./figures/module_signature_heatmap.pdf"),
     width = 50,height = 8)
draw(auc_signature)
dev.off()

saveRDS(ss1, "./data/rds/Step8_var4000.rds")

#####
#venndiagram
# Libraries
library(eulerr)
set.seed(1)

top3_wgcna_list=list("35.0_avp.crhb"=extracted_genes[["35.0_avp.crhb"]],
                "35.1_avp"=extracted_genes[["35.1_avp"]],
                "45_sst1.1"=extracted_genes[["45_sst1.1"]],
                #"41_nppal"=extracted_genes[["41_nppal"]],
                "27_Neuro"=extracted_genes[["27_Neuro"]]
                )

#a=euler(top3_wgcna_list,asp=1, shape ="circle" #ellipse circle )
#a.venn =  plot(a,quantities = T, fill = viridis(5), alpha =0.3,labels = "", lty =0)
#plot(a.venn)

a=ggvenn(top3_wgcna_list, show_elements = F, stroke_size = 0.5,stroke_alpha = 0.5,
         label_sep = "\n", text_size = 5,fill_alpha = 0.1,
         fill_color =  c("#4401541A","#21908D1A","#228B2233","#FFA50080" )
)
pdf(paste0("./figures/top3_WGCNA_venn.pdf"),
    width = 5,height = 5)

print(a)
dev.off() 


#GO enrichment
#####

#WGCNA module GO
# run_merged_WGCNA_GO(): GO enrichment over a merged/combined gene list across libraries
run_merged_WGCNA_GO = function(gene_list){
  library(xlsx)
  ## Gene Ontology Analysis
  #For this we will use a wrapper function of the GO enrichment analysis tools provided by the package `topGO`. First the setup.
  #{r danio_GOs_setup, echo = FALSE, warning = FALSE}
  #gene universe
  gene_universe <- rownames(ss1@assays$RNA$counts)
  
  #danio GOdb
  library(org.Dr.eg.db)
  allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Dr.eg.db", ID="symbol")
  # gene-GO mappings
  danio_id_GO <- allGO2genes
  
  #list of genes of interest
  danio_wg_list <-gene_list
  #And now for the GO analysis:
  
  #{r danio_GOs, message = FALSE}
  # GO term analysis wrapper
  danio_wg_GO_all<- 
    getGOs(
      danio_wg_list,
      gene_universe = gene_universe,
      gene2GO = danio_id_GO
    )
  
  # gene-GO mappings
  
  #brain background
  #ref_cell: markerset1_scheir_lab (Shafer et al., 2022, Nat Ecol Evol.)
  PO_scRNA_seq=readRDS("Z:/MinK/NGS_sequencing/sc_multiome_NPO_Project_11254/old_analysis/PO_scRNA_seq_v2.rds")
  expressed_genes_HYP <- rownames(PO_scRNA_seq)[Matrix::rowSums(PO_scRNA_seq@assays$RNA$counts > 0) > 0]
  
  # GO term analysis wrapper
  danio_wg_GO_all_rosetaall_bg_elim <- 
    getGOs(
      danio_wg_list,
      gene_universe = expressed_genes_HYP,
      gene2GO = danio_id_GO,
      alg = "elim"
    )
  
  pdf(paste0("./figures/danio_merged_wgcna_module_exploration_GOs_DIFFNORM_universe_rosettaBG_elim.pdf"))
  for (i in danio_wg_GO_all_rosetaall_bg_elim$GOplot) {print(i)}
  dev.off()
  
  # GO term analysis wrapper
  danio_wg_GO_all_rosetaall_bg <- 
    getGOs(
      danio_wg_list,
      gene_universe = expressed_genes_HYP,
      gene2GO = danio_id_GO
    )
  
  pdf(paste0("./figures/danio_merged_wgcna_module_exploration_GOs_DIFFNORM_universe_rosettaBG.pdf"))
  for (i in danio_wg_GO_all_rosetaall_bg$GOplot) {print(i)}
  dev.off()
  
  # GO term analysis wrapper
  danio_wg_GO_all_danioidmodule_bg <- 
    getGOs(
      danio_wg_list,
      gene_universe = danio_id_module$id,
      gene2GO = danio_id_GO
    )
  
  pdf(paste0("./figures/danio_merged_wgcna_module_exploration_GOs_DIFFNORM_universe_wgcnaBG.pdf"))
  for (i in danio_wg_GO_all_danioidmodule_bg$GOplot) {print(i)}
  dev.off()
  
  # GO term analysis wrapper
  danio_wg_GO_all_danioidmodule_bg_elim <- 
    getGOs(
      danio_wg_list,
      gene_universe = danio_id_module$id,
      gene2GO = danio_id_GO,
      alg = "elim"
    )
  
  pdf(paste0("./figures/danio_merged_wgcna_module_exploration_GOs_DIFFNORM_universe_wgcnaBG_elim.pdf"))
  for (i in danio_wg_GO_all_danioidmodule_bg_elim$GOplot) {print(i)}
  dev.off()
  
  
  
  #Here we show a couple of GO term analysis for different modules:
  
  #{r, fig.width = 6, fig.height = 6}
  danio_wg_GO_all$GOplot$s01
  
  pdf(paste0("./figures/danio_merged_wgcna_module_exploration_GOs_DIFFNORM_universe_allgenes.pdf"))
  for (i in danio_wg_GO_all$GOplot) {print(i)}
  dev.off()
  
  
  #####
  
  
  ## Saving genome region files for motif enrichment


## Saving everything

#{r save_GO_terms}
danio_wg_GO_table <- ldply(danio_wg_GO_all[[1]], .id="module")
write.xlsx(
  danio_wg_GO_table,
  file = paste0("./outputs/danio_merged_wgcna_GOterms.xlsx"),
  sheetName = "schmidtea_wgcna_GOterms",
  col.names = TRUE, row.names = FALSE, showNA = TRUE
)

pdf(paste0("./figures/danio_merged_wgcna_GOs_barplots.pdf"))
for (i in danio_wg_GO_all$GOplot) {print(i)}
dev.off()

danio_merged_wgcna_moduleinfo <- merge(danio_id_module, expressed_genes_HYP ,by = 1,)
write.table(
  danio_merged_wgcna_moduleinfo[danio_merged_wgcna_moduleinfo$gene_type == "hconf",],
  file = paste0("./outputs/danio_merged_wgcna_modules_info.tsv"),
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

}

run_merged_WGCNA_GO(extracted_genes)

####
#heatmap 


#####
#ATAC
DefaultAssay(ss1)="ATAC_macs3"

#extracted_genes=load("./data/rda/wgcna_genes.rda")
dir.create("./figures/ATAC/WGCNA")

#35.1_avp
avp_wgcna=extracted_genes[["35.1_avp"]]

for (i in avp_wgcna) {
  gene_cordi=LookupGeneCoords(ss1,i)
  if(is.null(gene_cordi) == FALSE){
    hits <- findOverlaps(ranges.show, gene_cordi)
    
    # Extract elements that fall within the range
    selected_elements <- ranges.show[queryHits(hits)]
    subgroup2show="35.1_avp"
    a=CoveragePlot(
      object = ss1,split.by = "orig.ident",#peaks.group.by = "orig.ident",  
      group.by = "merged_sub.anno_type",  
      region = gene_cordi,
      features = i,
      region.highlight = selected_elements,
      links = T,
      #heights = c(18,1,1,1),
      expression.assay = "RNA",
      idents = subgroup2show,
      extend.upstream = 2000,
      extend.downstream = 2000
    )
    
    
    n=str_split(str_split(string = a[[1]][[1]][["labels"]][["y"]],pattern = "- ")[[1]][2],"\\)")[[1]][1]
    if (as.integer(n) > 4) {
      tiff(paste0("./figures/ATAC/WGCNA/ATAC_avp_",i,".tiff"),
           width = 15,height = 20,units = "cm", res = 300,compression = "lzw")
      print(a & scale_fill_manual(values = magma(5,alpha = 0.5)))
      dev.off()
      aa=paste0(i," has enough peak")
      print(aa)
    }else{
      aa=paste0(i," not enough peak")
      print(aa)
    }
  }
}

#35.0_avp.crhb
avp.crhb_wgcna=extracted_genes[["35.0_avp.crhb"]]

for (i in avp_wgcna) {
  gene_cordi=LookupGeneCoords(ss1,i)
  if(is.null(gene_cordi) == FALSE){
    hits <- findOverlaps(ranges.show, gene_cordi)
    
    # Extract elements that fall within the range
    selected_elements <- ranges.show[queryHits(hits)]
    subgroup2show="35.0_avp.crhb"
    a=CoveragePlot(
      object = ss1,split.by = "orig.ident",#peaks.group.by = "orig.ident",  
      group.by = "merged_sub.anno_type",  
      region = gene_cordi,
      features = i,
      region.highlight = selected_elements,
      links = T,
      #heights = c(18,1,1,1),
      expression.assay = "RNA",
      idents = subgroup2show,
      extend.upstream = 2000,
      extend.downstream = 2000
    )
    
    
    n=str_split(str_split(string = a[[1]][[1]][["labels"]][["y"]],pattern = "- ")[[1]][2],"\\)")[[1]][1]
    if (as.integer(n) > 4) {
      tiff(paste0("./figures/ATAC/WGCNA/ATAC_avp.crhb_",i,".tiff"),
           width = 15,height = 20,units = "cm", res = 300,compression = "lzw")
      print(a & scale_fill_manual(values = magma(5,alpha = 0.5)))
      dev.off()
      aa=paste0(i," has enough peak")
      print(aa)
    }else{
      aa=paste0(i," not enough peak")
      print(aa)
    }
  }
}

#sst
sst_wgcna=extracted_genes[["45_sst1.1"]]

for (i in sst_wgcna) {
  gene_cordi=LookupGeneCoords(ss1,i)
  if(is.null(gene_cordi) == FALSE){
    hits <- findOverlaps(ranges.show, gene_cordi)
    
    # Extract elements that fall within the range
    selected_elements <- ranges.show[queryHits(hits)]
    subgroup2show="45_sst1.1"
    
    a=CoveragePlot(
      object = ss1,split.by = "orig.ident",#peaks.group.by = "orig.ident",  
      group.by = "merged_sub.anno_type",  
      region = gene_cordi,
      features = i,
      region.highlight = selected_elements,
      links = T,
      #heights = c(18,1,1,1),
      expression.assay = "RNA",
      idents = subgroup2show,
      extend.upstream = 2000,
      extend.downstream = 2000
    )
    
    
    n=str_split(str_split(string = a[[1]][[1]][["labels"]][["y"]],pattern = "- ")[[1]][2],"\\)")[[1]][1]
    if (as.integer(n) > 4) {
      tiff(paste0("./figures/ATAC/WGCNA/ATAC_sst_",i,".tiff"),
           width = 15,height = 20,units = "cm", res = 300,compression = "lzw")
      print(a & scale_fill_manual(values = magma(5,alpha = 0.5)))
      dev.off()
      aa=paste0(i," has enough peak")
      print(aa)
    }else{
      aa=paste0(i," not enough peak")
      print(aa)
    }
  }
}


#avp, avp.crhb. sst
#####

cell_group=c("35.0_avp.crhb", "35.1_avp","45_sst1.1")


for (k in 1:3) {
  
  wgcnas=extracted_genes[[(cell_group[k])]]

  for (i in wgcnas) {
    gene_cordi=LookupGeneCoords(ss1,i)
    if(is.null(gene_cordi) == FALSE){
      hits <- findOverlaps(ranges.show, gene_cordi)
      
      # Extract elements that fall within the range
      selected_elements <- ranges.show[queryHits(hits)]
      
      a=CoveragePlot(
        object = ss1,
        #split.by = "orig.ident",#peaks.group.by = "orig.ident",  
        group.by = "merged_sub.anno_type",  
        region = gene_cordi,
        features = i,
        region.highlight = selected_elements,
        links = T,
        #heights = c(18,1,1,1),
        expression.assay = "RNA",
        idents = cell_group,
        extend.upstream = 2000,
        extend.downstream = 2000
      )
      
      
      n=str_split(str_split(string = a[[1]][[1]][["labels"]][["y"]],pattern = "- ")[[1]][2],"\\)")[[1]][1]
      if (as.integer(n) > 4) {
        tiff(paste0("./figures/ATAC/WGCNA/ATAC_np3_",cell_group[k],"_",i,".tiff"),
             width = 15,height = 15,units = "cm", res = 300,compression = "lzw")
        print(a & scale_fill_manual(values = magma(5,alpha = 0.5)))
        dev.off()
        aa=paste0(i," has enough peak")
        print(aa)
      }else{
        aa=paste0(i," not enough peak")
        print(aa)
      }
    }
  }
  
}

##### signac motif : TF-motif enrichment in the accessible regions of module genes
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Drerio.UCSC.danRer11)
library(patchwork)


pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
DefaultAssay(ss1)="ATAC_macs3"
ss1 <- AddMotifs(
  object = ss1,
  genome = BSgenome.Drerio.UCSC.danRer11,
  pfm = pfm
)

saveRDS(ss1,paste0("./data/rds/Step8_motif_var",r.variable,".rds"))
'
ss1=readRDS(paste0("./data/rds/Step8_motif_var",r.variable,".rds"))
library(rtracklayer)
export.bed(ss1@assays$ATAC_macs3@annotation,
           con="./outputs/peak_annotation_granges.bed")
           

export.bed(granges(ss1[["ATAC_macs3"]]), "./outputs/ATAC_macs3_peaks.bed")

#ptp4a2b

a=FindMotifs(
    ss1,
    features="chr19-35455566-35456319",
    background = 50000,
    assay = "ATAC_macs3",
    verbose = TRUE,
    p.adjust.method = "BH"
    
)


aa=dplyr::filter(a,observed >0)


aaa=MotifPlot(
  object = ss1,
  motifs = head(rownames(enriched.motifs),20)
)



# find peaks open in cells

Idents(ss1)=ss1$merged_sub.anno_type_ori
open.peaks <- AccessiblePeaks(ss1, 
                              idents = c("45_sst1.1:cont_post", "45_sst1.1:cont_pre"),
                              assay = "ATAC_macs3")

# match the overall GC content in the peak set
meta.feature <- GetAssayData(ss1, assay = "ATAC_macs3", layer = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  n = 50000
)
'

#bPAC
da_peaks <- FindMarkers(
  object = ss1,
  #reduction = "wnn.umap",
  assay = "ATAC_macs3",
  group.by = "merged_sub.anno_type_ori",
  ident.1 = '45_sst1.1:bPAC_post',
  ident.2 = '45_sst1.1:bPAC_pre',#'45_sst1.1:cont_post',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_ATAC_macs3'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.025 & da_peaks$pct.1 > 0.2, ])  
top.da.peak.loc <- da_peaks[da_peaks$p_val < 0.025 & da_peaks$pct.1 > 0.2, ]
top.peak.links=ss1@assays$ATAC_macs3@links[which(ss1@assays$ATAC_macs3@links$peak %in% top.da.peak )]
write.csv(as.data.table(top.peak.links),"./outputs/ATAC_res/top.peak.links_bPAC.csv")

sst.top.peak=unique(top.peak.links$gene)

#control
da_peaks_cont <- FindMarkers(
  object = ss1,
  #reduction = "wnn.umap",
  assay = "ATAC_macs3",
  group.by = "merged_sub.anno_type_ori",
  ident.1 = '45_sst1.1:cont_post',
  ident.2 = '45_sst1.1:cont_pre',#'45_sst1.1:cont_post',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_ATAC_macs3'
)

# get top differentially accessible peaks
top.da.peak_cont <- rownames(da_peaks_cont[da_peaks_cont$p_val < 0.025 & da_peaks_cont$pct.1 > 0.2, ])  
top.da.peak_cont.loc <- da_peaks_cont[da_peaks_cont$p_val < 0.025 & da_peaks_cont$pct.1 > 0.2, ]
top.da.peak_cont.links=ss1@assays$ATAC_macs3@links[which(ss1@assays$ATAC_macs3@links$peak %in% top.da.peak_cont )]
write.csv(as.data.table(top.da.peak_cont.links),"./outputs/ATAC_res/top.peak.links_cont.csv")

sst.top.peak_cont=unique(top.da.peak_cont.links$gene)


intersect(sst.top.peak_cont,sst.top.peak)


deg45_bPAC=read.csv("./outputs/DEGs/DEGs_postLD_preLD_45.sst1.1.bPAC.post_45.sst1.1.bPAC.pre.csv",row.names = 1)
deg45_bPAC=dplyr::filter(deg45_bPAC, padj <0.05)

deg45_post=read.csv("./outputs/DEGs/DEGs_bPAC_cont_45.sst1.1.bPAC.post_45.sst1.1.cont.post.csv",row.names = 1)
deg45_post=dplyr::filter(deg45_post, pvalue <0.005)

intersect(sst.top.peak, rownames(deg45_bPAC))
intersect(sst.top.peak, rownames(deg45_post)) 
#bPAC
enriched.motifs <- FindMotifs(
  object = ss1,
  features = top.da.peak,
  background = 40000,
  assay = "ATAC_macs3",
  verbose = TRUE,
  p.adjust.method = "BH"
  
)

#enriched.motifs=dplyr::filter(enriched.motifs,p.adjust <0.01)

a=MotifPlot(
  object = ss1,
  motifs = head(rownames(enriched.motifs),20)
)

#control
enriched.motifs_cont <- FindMotifs(
  object = ss1,
  features = top.da.peak_cont,
  background = 40000,
  assay = "ATAC_macs3",
  verbose = TRUE,
  p.adjust.method = "BH"
  
)

#enriched.motifs_cont=dplyr::filter(enriched.motifs_cont,p.adjust <0.01)

aa=MotifPlot(
  object = ss1,
  motifs = head(rownames(enriched.motifs_cont),20)
)
pdf(paste0("./figures/ATAC/top_peaks/45_sst_enriched_motifs_cont_bPAC.pdf"),
    width = 10,height =10)
plot(aa/a)
dev.off()
library(ggrepel)

DMR_enriched_TF=full_join(enriched.motifs_cont,enriched.motifs,by="motif.name")
DMR_enriched_TF_tb= DMR_enriched_TF[,c(8,6,15,9,17)]
DMR_enriched_TF_tb$fold.enrichment.x[is.na(DMR_enriched_TF$fold.enrichment.x)]=0
DMR_enriched_TF_tb$fold.enrichment.y[is.na(DMR_enriched_TF$fold.enrichment.y)]=0
DMR_enriched_TF_tb$p.adjust.x[is.na(DMR_enriched_TF$p.adjust.x)]=1
DMR_enriched_TF_tb$p.adjust.y[is.na(DMR_enriched_TF$p.adjust.y)]=1

DMR_enriched_TF_tb[,"delta_FC"]=DMR_enriched_TF_tb$fold.enrichment.y-DMR_enriched_TF_tb$fold.enrichment.x
labels_pt=DMR_enriched_TF_tb$motif.name
labels_pt[which(DMR_enriched_TF_tb$p.adjust.x >0.0005 & DMR_enriched_TF_tb$p.adjust.y >0.0005)]=""
labels_pt[which(abs(DMR_enriched_TF_tb$delta_FC)<quantile(abs(DMR_enriched_TF_tb$delta_FC))[4])]=""
DMR_enriched_TF_tb[,"labels"]=labels_pt

DMR_enriched_TF_tb[,"sum_FC"]=DMR_enriched_TF_tb$fold.enrichment.y+DMR_enriched_TF_tb$fold.enrichment.x
labels_pt2=DMR_enriched_TF_tb$motif.name
labels_pt2[which(DMR_enriched_TF_tb$p.adjust.x >0.0005 | DMR_enriched_TF_tb$p.adjust.y >0.0005)]=""
labels_pt2[which(abs(DMR_enriched_TF_tb$sum_FC)<quantile(abs(DMR_enriched_TF_tb$sum_FC))[4])]=""
DMR_enriched_TF_tb[,"labels2"]=labels_pt2

DMR_enriched_TF_tb[,4]=-log10(DMR_enriched_TF_tb[,4])
DMR_enriched_TF_tb[,5]=-log10(DMR_enriched_TF_tb[,5])


p <- ggplot(DMR_enriched_TF_tb, aes(p.adjust.y,p.adjust.x ,label=labels, color=delta_FC))

pdf(paste0("./figures/ATAC/top_peaks/45_sst_enriched_motif_comp.pdf"),
    width = 6.5,height =6)
p + geom_point()+
  geom_abline(slope = 1, intercept = 0,linetype=3)+
  geom_hline(yintercept = -log10(0.0005),linetype=3)+
  geom_vline(xintercept = -log10(0.0005),linetype=3)+
  geom_text_repel()+
  ylim(c(0,7.5))+xlim(c(0,7.5))+
  theme_classic()+
  scale_colour_gradientn(
    colours = c("blue", "dodgerblue", "lightgray", "orange", "red"),
    rescaler = ~ scales::rescale_mid(.x, mid = 0)
  )
dev.off()

p <- ggplot(DMR_enriched_TF_tb, aes(p.adjust.y,p.adjust.x ,label=labels2, color=scale(sum_FC)))

pdf(paste0("./figures/ATAC/top_peaks/45_sst_enriched_motif_comp2.pdf"),
    width = 6.5,height =6)
p + geom_point()+
  geom_abline(slope = 1, intercept = 0,linetype=3)+
  geom_hline(yintercept = -log10(0.0005),linetype=3)+
  geom_vline(xintercept = -log10(0.0005),linetype=3)+
  geom_text_repel()+
  ylim(c(0,7.5))+xlim(c(0,7.5))+
  theme_classic()+
  scale_colour_gradientn(
    colours = c("blue", "dodgerblue", "lightgray", "orange", "red"),
    rescaler = ~ scales::rescale_mid(.x, mid = 0)
  )
dev.off()



DMR_enriched_TF_tb[,"ratio"]= DMR_enriched_TF_tb$p.adjust.y/DMR_enriched_TF_tb$p.adjust.x
DMR_enriched_TF_bPAc=DMR_enriched_TF_tb[c(which(DMR_enriched_TF_tb$ratio >= 1.5)),]


save(da_peaks,da_peaks_cont,enriched.motifs,enriched.motifs_cont,
     file = paste0("./data/rda/enrichedmotifs.rda")
)



#####

#ptp4a2b

ptp4a2b_peak1.motif=FindMotifs(
  ss1,
  features="chr19-35455566-35456319",
  background = 50000,
  assay = "ATAC_macs3",
  verbose = TRUE,
  p.adjust.method = "BH"
  
)

ft.ptp4a2b_peak1.motifa= ptp4a2b_peak1.motif %>% dplyr::filter(observed >0) %>% dplyr::filter(pvalue < 0.1)



#DEGs in 45
deg45_bPAC=read.csv("./outputs/DEGs/DEGs_postLD_preLD_45.sst1.1.bPAC.post_45.sst1.1.bPAC.pre.csv",row.names = 1)
deg45_bPAC=dplyr::filter(deg45_bPAC, pvalue <0.005)
TFs=na.omit(data.frame("gene"=rownames(ss1@assays$RNA) ,"tf"=ss1@assays$RNA@meta.features$TF))


MotifPlot(
  object = ss1,
  motifs = head(rownames(ft.ptp4a2b_peak1.motifa),12)
)



#DEGs

topLD_cluster=c("35.0.avp.crhb","35.1.avp","45.sst1.1" )
DEG_lists=list()
for (i in topLD_cluster) {
  deg_bPAC=read.csv(paste0("./outputs/DEGs/DEGs_postLD_preLD_",i,".bPAC.post_",i,".bPAC.pre.csv"),row.names = 1)
  deg_bPAC=dplyr::filter(deg_bPAC, padj <0.1)
  deg_cont=read.csv(paste0("./outputs/DEGs/DEGs_postLD_preLD_",i,".cont.post_",i,".cont.pre.csv"),row.names = 1)
  deg_cont=dplyr::filter(deg_cont, padj <0.1)
  deg_pre=read.csv(paste0("./outputs/DEGs/DEGs_bPAC_cont_",i,".bPAC.pre_",i,".cont.pre.csv"),row.names = 1)
  deg_pre=dplyr::filter(deg_pre, padj <0.1)
  deg_post=read.csv(paste0("./outputs/DEGs/DEGs_bPAC_cont_",i,".bPAC.post_",i,".cont.post.csv"),row.names = 1)
  deg_post=dplyr::filter(deg_post, padj <0.1)
  
  DEG_lists[paste0("deg_bPAC_",i)]= list(deg_bPAC)
  DEG_lists[paste0("deg_cont_",i)]= list(deg_cont)
  DEG_lists[paste0("deg_pre_",i)]= list(deg_pre)
  DEG_lists[paste0("deg_post_",i)]= list(deg_post)
}


#expressed
all_marks=read.csv("./outputs/cell_type/all.sub.markers_var4000_merged_sub.anno.csv",row.names = 1)
c45_mark=all_marks %>% dplyr::filter(cluster == "45_sst1.1")%>% dplyr::filter(p_val_adj <0.05)
TFs=na.omit(data.frame("gene"=rownames(ss1@assays$RNA) ,"tf"=ss1@assays$RNA@meta.features$TF))


en.motif= unique(c(sapply(str_split(string = tolower(unique(ft.ptp4a2b_peak1.motifa$motif.name)),
                                    pattern = "\\(|\\:"), `[`, 1),na.omit(sapply(str_split(string = tolower(unique(ft.ptp4a2b_peak1.motifa$motif.name)), pattern = "\\(|\\:"), `[`, 3)) ))
overlapped= list()

for (i in en.motif) {
  overlapped[[i]]=c45_mark$gene[grepl(pattern = i, x = c45_mark$gene)]
}

#overlapped
exp_enriched_motifs=toupper(c("creb3l1","erf","rora","fos","fosl2","fosab","fosb","jun","jund"))
exp_motif.ptp4a2b= c()
for (i in exp_enriched_motifs) {
  exp_motif.ptp4a2b=c(exp_motif.ptp4a2b,ft.ptp4a2b_peak1.motifa$motif.name[grep(i,ft.ptp4a2b_peak1.motifa$motif.name)])
}
exp_motif.ptp4a2b=unique(exp_motif.ptp4a2b)

#exp_enriched_motifs=c("NR3C2","Ar","CREM","Foxo1","FOSL2::JUN(var.2)","FOSL2::JUND(var.2)","FOSL2::JUNB(var.2)","KLF9","SP4")

pdf(paste0("./figures/ATAC/top_peaks/45_sst_exp_ptp42b_motif.pdf"),
    width = 15,height = 6)
MotifPlot(
  object = ss1,
  motifs = head(rownames(ft.ptp4a2b_peak1.motifa[which(ft.ptp4a2b_peak1.motifa$motif.name %in% exp_motif.ptp4a2b),]),24), ncol=4
)
dev.off()



FeatureDimPlot(
  srt = ss1, 
  features = c("nr3c2",  "ar", "foxo1"),
  assay = "SCT",
  seed = 0, compare_features = F, label =F, label_repel = T,label_insitu = TRUE, 
  add_density = F,palette = "GdRd",bg_cutoff = 0.5, 
  reduction = "wnn.umap", 
  theme_use = "theme_blank",
  theme_args=theme(strip.text = element_text(size = 20,face = c("bold.italic")))
)

dir.create("./figures/ATAC/top_peaks")
for (i in sst.top.peak) {
  
  subgroup2show=c("35.0_avp.crhb","35.1_avp","45_sst1.1" )
  
  a=CoveragePlot(
    object = ss1,split.by = "orig.ident",peaks.group.by = "merged_sub.anno_type",  
    group.by = "merged_sub.anno_type",  
    region = i,
    features = i,
    #region.highlight = selected_elements,
    links = T,
    #heights = c(18,1,1,1),
    expression.assay = "RNA",
    idents = subgroup2show,
    extend.upstream = 2000,
    extend.downstream = 2000
  )
  
  
  n=str_split(str_split(string = a[[1]][[1]][["labels"]][["y"]],pattern = "- ")[[1]][2],"\\)")[[1]][1]
  if (as.integer(n) > 4) {
    pdf(paste0("./figures/ATAC/top_peaks/ATAC_sst_",i,".pdf"),
        width = 15,height = 20)
    print(a & scale_fill_manual(values = magma(13,alpha = 0.5)))
    dev.off()
    aa=paste0(i," has enough peak")
    print(aa)
  }else{
    aa=paste0(i," not enough peak")
    print(aa)
  }
}

#s.fig11

#DEG_lists
DEG45_sst.top.peak=intersect(sst.top.peak, rownames(DEG_lists$deg_bPAC_45.sst1.1))

sst_deg_top_cover_plot=list()
for (i in DEG45_sst.top.peak) {
  
  subgroup2show=c("45_sst1.1" )
  
  a=CoveragePlot(
    object = ss1,split.by = "orig.ident",#peaks.group.by = "merged_sub.anno_type",  
    group.by = "merged_sub.anno_type",  
    region = i,
    features = i,
    #region.highlight = selected_elements,
    links = T,
    #heights = c(4,2,2,2),
    expression.assay = "RNA",
    idents = subgroup2show,
    extend.upstream = 2000,
    extend.downstream = 2000
  )
  
  
  n=str_split(str_split(string = a[[1]][[1]][["labels"]][["y"]],pattern = "- ")[[1]][2],"\\)")[[1]][1]
  if (as.integer(n) > 4) {
    aa=paste0(i," has enough peak")
    print(aa)
    sst_deg_top_cover_plot[i]=a & scale_fill_manual(values = magma(5,alpha = 0.5))
  }
}


pdf(paste0("./figures/ATAC/top_peaks/ATAC_sst_45_DEG.pdf"),
    width = 30,height = 15)
wrap_plots(sst_deg_top_cover_plot, ncol = 4)
dev.off()

