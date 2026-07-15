# =============================================================================
# 04_1_ATAC_annotation_11_2025v3.R
# -----------------------------------------------------------------------------
# Purpose : Link MACS3 ATAC peaks to genes via peak-to-expression correlation
#           (Signac LinkPeaks), producing candidate cis-regulatory elements for
#           each gene. Adds sequence stats needed for the linkage model.
# Inputs  : Seurat object `ss1` (must contain an "ATAC_macs3" peak assay + RNA).
#           AnnotationHub EnsDb (Danio rerio, Ensembl 111); danRer11 genome.
# Output  : ./data/rds/Step4_var<r.variable>.rds  (object with peak-gene links)
# =============================================================================

## Load libraries
library(plyr)
library(dplyr)
library(ggplot2)
library(Seurat)
library(colorspace)
library(viridis)
library(SCP)
#set the data and annotations
set.seed(1234)
options(future.globals.maxSize = 120000 * 1024^2)   # raise memory cap for large objects

r.variable=4000       # number of variable features (encoded in RDS names)
vs="06_2025_v4"       # version tag
umap="wnn.umap"       # UMAP embedding to use for plots

#set working dir (skip if already inside the version folder)
if(getwd()!=paste0("D:/projects/scMultiome_oxt/",vs)){
  setwd(paste0("./",vs))
}

#functions adopted and modified from Alberto Perez-Posada @apposada
source("./ext_code/r_code/functions/sourcefolder.R")
sourceFolder("./ext_code/r_code/functions/")

#other temp function: negated %in%
'%!in%' <- function(x,y)!('%in%'(x,y))
#load data
#ss1=readRDS(paste0("./data/rds/step4_var",r.variable,".rds"))

# ATAC_link ---------------------------------------------------------------
# Build zebrafish gene annotation (UCSC/danRer11 style) for peak-gene linkage.

library(Signac)
library(Seurat)
library(BSgenome.Drerio.UCSC.danRer11)
library(AnnotationHub)

#library(biomaRt)
#library(org.Dr.eg.db)

ah = AnnotationHub()
ensdbs <- query(ah, c("Danio rerio"))
ensdb_id <- ensdbs$ah_id[grep(paste0(" 111 EnsDb"), ensdbs$title)]   # Ensembl release 111
ensdb <- ensdbs[[ensdb_id]]
annotations <- GetGRangesFromEnsDb(ensdb = ensdb, standard.chromosomes =T)
seqlevelsStyle(annotations) = "UCSC"
genome(annotations) = "danRer11"

#load seurat obj

DefaultAssay(ss1) <- "ATAC_macs3"

# RegionStats adds GC content / sequence stats per peak (required by LinkPeaks)
ss1 <- RegionStats(ss1, genome = BSgenome.Drerio.UCSC.danRer11)

# link peaks to genes: correlate peak accessibility with gene expression across cells
ss1 <- LinkPeaks(
  object = ss1,
  peak.assay = "ATAC_macs3",
  expression.assay = "RNA",

)

# save Seurat objects -----------------------------------------------------

saveRDS(ss1,paste0("./data/rds/Step4_var",r.variable,".rds"))




