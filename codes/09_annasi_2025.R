# =============================================================================
# 09_annasi_2025.R
# -----------------------------------------------------------------------------
# Purpose : Export pseudobulk RNA + ATAC matrices and contrasts for scANANSE /
#           ANANSE gene-regulatory-network inference (via AnanseSeurat). Splits
#           each cell group into two pseudo-replicates (R1/R2), then writes the
#           CPM, ATAC, config, and DEG files ANANSE needs per cell type/condition.
# Inputs  : ./data/rds/Step7_var<r.variable>.rds
#           ./outputs/cell_type/cell_type_table_m_modified.csv
# Outputs : ./scANANSE/analysis_11_2025/  (scANANSE input files + contrasts)
# Note    : "annasi"/"ANANSE" = ANalysis Algorithm for Networks Specified by
#           Enhancers; downstream network building is in the ananse_graph_* files.
# =============================================================================

library(Seurat)
library(SeuratDisk)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(AnanseSeurat)      # bridges Seurat objects to scANANSE exports
library(SeuratDisk)

#set the data and annotations
set.seed(1234)
options(future.globals.maxSize = 120000 * 1024^2)

r.variable=4000
vs="06_2025_v4"
umap="wnn.umap"

#set working dir (skip if already inside the version folder)
if(getwd()!=paste0("D:/projects/scMultiome_oxt/",vs)){
  setwd(paste0("./",vs))
}

#functions adopted and modified from Alberto Perez-Posada @apposada
source("./ext_code/r_code/functions/sourcefolder.R")
sourceFolder("./ext_code/r_code/functions/")

#other temp function: negated %in%
'%!in%' <- function(x,y)!('%in%'(x,y))

## Load Data

ss1=readRDS(paste0("./data/rds/Step7_var",r.variable,".rds"))
libs=unique(ss1$orig.ident)
type_table_m.v2= read.csv(paste0("./outputs/cell_type/cell_type_table_m_modified.csv"), header = T, row.names = 1)
type_table_m.v3= type_table_m.v2[- c(grep(pattern = "UnD", type_table_m.v2$numbered)),]  # drop undefined types

#wt only

## Getting cluster information for pseudobulk

#scRNA-seq:
DefaultAssay(ss1)<-"RNA"

#pseudobulk for condition x clusters
#ss1$merged_sub.anno_type_ori=paste0(ss1$merged_sub.anno_type,".",ss1$orig.ident)

Idents(ss1)=ss1$modifID_cc                 # identity = cell type x condition label

#ss1$modifID_cc=gsub("\\.","-",as.character(ss1$modifID_cc))

# pseudo_split ------------------------------------------------------------
# ANANSE needs replicates. Randomly split each group's cells into two halves
# (R1/R2) to create two pseudo-replicates per cell type/condition.

df <- data.frame(ID = names(ss1$modifID_cc), group = ss1$modifID_cc)

# Split each group randomly into two, preserving IDs
set.seed(123)
split_groups <- split(df, df$group)

group1 <- c()
group2 <- c()

for (grp in split_groups) {
  n <- nrow(grp)
  idx <- sample(1:n, floor(n / 2))        # random half -> replicate 1
  group1 <- c(group1, paste0(grp$ID[idx]))
  group2 <- c(group2, paste0(grp$ID[-idx]))  # remaining half -> replicate 2
}

pseudo_anno=ss1$modifID_cc
pseudo_anno[group1]=paste0(pseudo_anno[group1],".R1")
pseudo_anno[group2]=paste0(pseudo_anno[group2],".R2")

ss1$pseudo_anno=pseudo_anno               # per-cell replicate-tagged label
Idents(ss1)=ss1$pseudo_anno

#### Build pre- vs post-condition contrasts for each cell type
post_cyte=sort(unique(ss_idents)[grep("post",unique(ss_idents))])
pre_cyte=sort(unique(ss_idents)[grep("pre",unique(ss_idents))])


list_contrasts <-
  as.list(as.data.frame(t(data.frame(
    "ctype",post_cyte,pre_cyte
  ))))


names(list_contrasts) <- sapply(list_contrasts, function(x){paste(x[c(2,3)], collapse = "_")})

list_contrasts

# Reformat contrasts to ANANSE's expected "group1_group2" naming (dots not underscores)
add_contrasts=c()
for (i in names(list_contrasts)) {
  add_contrasts=c(add_contrasts,paste0(gsub("_",replacement = ".",list_contrasts[[i]][2]),
                                       "_",
                                       gsub("_",replacement = ".",list_contrasts[[i]][3])))
}



# exporting ---------------------------------------------------------------
# Write the scANANSE input bundle: CPM expression, ATAC peaks, config, and DEGs.

#coi= sort(unique(ss1$modifID_cc[!grepl("sub",ss1$modifID_cc)]))
dir.create("./scANANSE/")
outdir="./scANANSE/analysis_11_2025"
dir.create(outdir)
# Pseudobulk RNA counts per pseudo-replicate -> CPM table
export_CPM_scANANSE(
  ss1,
  min_cells = 4,
  output_dir = outdir,
  cluster_id = "pseudo_anno",
  RNA_count_assay = 'RNA'
)

# Pseudobulk ATAC (MACS3 peaks) per pseudo-replicate
export_ATAC_scANANSE(
  ss1,
  min_cells = 4,
  output_dir = outdir,
  cluster_id = "pseudo_anno",
  ATAC_peak_assay = 'ATAC_macs3'
)

# Specify additional contrasts:


add.contrasts <- add_contrasts


# Write the scANANSE config (genome, contrasts) at the cell-type level
config_scANANSE(
  ss1, genome = "danRer11",
  min_cells = 4,
  output_dir = outdir,
  cluster_id = "modifID_cc",
  additional_contrasts =add.contrasts
)


# Differential expression per contrast, fed into ANANSE influence scoring
DEGS_scANANSE(
  ss1,
  genome = "danRer11",
  min_cells = 4,
  output_dir = outdir,
  cluster_id = "modifID_cc",
  additional_contrasts = add.contrasts
)

