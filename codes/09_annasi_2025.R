library(Seurat)
library(SeuratDisk)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(AnanseSeurat)
library(SeuratDisk)

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

## Load Data

ss1=readRDS(paste0("./data/rds/Step7_var",r.variable,".rds"))
libs=unique(ss1$orig.ident)
type_table_m.v2= read.csv(paste0("./outputs/cell_type/cell_type_table_m_modified.csv"), header = T, row.names = 1)
type_table_m.v3= type_table_m.v2[- c(grep(pattern = "UnD", type_table_m.v2$numbered)),]

#wt only

## Getting cluster information for pseudobulk

#scRNA-seq:
DefaultAssay(ss1)<-"RNA"

#pseudobulk for condition x clusters
#ss1$merged_sub.anno_type_ori=paste0(ss1$merged_sub.anno_type,".",ss1$orig.ident)

Idents(ss1)=ss1$modifID_cc

#ss1$modifID_cc=gsub("\\.","-",as.character(ss1$modifID_cc))

# pseudo_split ------------------------------------------------------------

df <- data.frame(ID = names(ss1$modifID_cc), group = ss1$modifID_cc)

# Split each group randomly into two, preserving IDs
set.seed(123)
split_groups <- split(df, df$group)

group1 <- c()
group2 <- c()

for (grp in split_groups) {
  n <- nrow(grp)
  idx <- sample(1:n, floor(n / 2))
  group1 <- c(group1, paste0(grp$ID[idx]))
  group2 <- c(group2, paste0(grp$ID[-idx]))
}

pseudo_anno=ss1$modifID_cc
pseudo_anno[group1]=paste0(pseudo_anno[group1],".R1")
pseudo_anno[group2]=paste0(pseudo_anno[group2],".R2")

ss1$pseudo_anno=pseudo_anno
Idents(ss1)=ss1$pseudo_anno

####
post_cyte=sort(unique(ss_idents)[grep("post",unique(ss_idents))])
pre_cyte=sort(unique(ss_idents)[grep("pre",unique(ss_idents))])


list_contrasts <-
  as.list(as.data.frame(t(data.frame(
    "ctype",post_cyte,pre_cyte
  ))))


names(list_contrasts) <- sapply(list_contrasts, function(x){paste(x[c(2,3)], collapse = "_")})

list_contrasts

add_contrasts=c()
for (i in names(list_contrasts)) {
  add_contrasts=c(add_contrasts,paste0(gsub("_",replacement = ".",list_contrasts[[i]][2]),
                                       "_",
                                       gsub("_",replacement = ".",list_contrasts[[i]][3])))
}



# exporting ---------------------------------------------------------------

#coi= sort(unique(ss1$modifID_cc[!grepl("sub",ss1$modifID_cc)]))
dir.create("./scANANSE/")
outdir="./scANANSE/analysis_11_2025"
dir.create(outdir)
export_CPM_scANANSE(
  ss1,
  min_cells = 4,
  output_dir = outdir,
  cluster_id = "pseudo_anno",
  RNA_count_assay = 'RNA'
)

export_ATAC_scANANSE(
  ss1,
  min_cells = 4,
  output_dir = outdir,
  cluster_id = "pseudo_anno",
  ATAC_peak_assay = 'ATAC_macs3'
)

# Specify additional contrasts:


add.contrasts <- add_contrasts


config_scANANSE(
  ss1, genome = "danRer11",
  min_cells = 4,
  output_dir = outdir,
  cluster_id = "modifID_cc",
  additional_contrasts =add.contrasts
)


DEGS_scANANSE(
  ss1,
  genome = "danRer11",
  min_cells = 4,
  output_dir = outdir,
  cluster_id = "modifID_cc",
  additional_contrasts = add.contrasts
)

