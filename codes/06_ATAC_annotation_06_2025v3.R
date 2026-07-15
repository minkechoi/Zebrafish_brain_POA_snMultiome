# =============================================================================
# 06_ATAC_annotation_06_2025v3.R
# -----------------------------------------------------------------------------
# Purpose : Peak-to-gene linkage (Signac LinkPeaks) for two objects:
#           (1) the full dataset `ss1` and (2) the neuropeptidergic subset `ss2`
#           (nps = neuropeptide-expressing neurons). Same procedure as 04_1 but
#           applied to both the whole atlas and the focused subset.
# Inputs  : `ss1` in memory; ./data/rds/step5_nps_var<r.variable>.rds for ss2.
#           AnnotationHub EnsDb (Danio rerio, Ensembl 111); danRer11 genome.
# Outputs : ./data/rds/Step6_var<r.variable>.rds       (full, with links)
#           ./data/rds/Step6_nps_var<r.variable>.rds   (nps subset, with links)
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
#load data
#ss1=readRDS(paste0("./data/rds/step4_var",r.variable,".rds"))

# ATAC_link ---------------------------------------------------------------
# Build zebrafish gene annotation (UCSC/danRer11) for peak-gene linkage.

library(Signac)
library(Seurat)
library(BSgenome.Drerio.UCSC.danRer11)
library(AnnotationHub)

#library(biomaRt)
#library(org.Dr.eg.db)

ah = AnnotationHub()
ensdbs <- query(ah, c("Danio rerio"))
ensdb_id <- ensdbs$ah_id[grep(paste0(" 111 EnsDb"), ensdbs$title)]
ensdb <- ensdbs[[ensdb_id]]
annotations <- GetGRangesFromEnsDb(ensdb = ensdb, standard.chromosomes =T)
seqlevelsStyle(annotations) = "UCSC"
genome(annotations) = "danRer11"

#load seurat obj

DefaultAssay(ss1) <- "ATAC_macs3"

# Add per-peak sequence stats (required by LinkPeaks) then link peaks to genes
ss1 <- RegionStats(ss1, genome = BSgenome.Drerio.UCSC.danRer11)

# link peaks to genes (full dataset)
ss1 <- LinkPeaks(
  object = ss1,
  peak.assay = "ATAC_macs3",
  expression.assay = "RNA",

)

# save Seurat objects -----------------------------------------------------

saveRDS(ss1,paste0("./data/rds/Step6_var",r.variable,".rds"))

##### for npc  (neuropeptidergic-neuron subset)

ss2= readRDS(paste0("./data/rds/step5_nps_var",r.variable,".rds"))
DefaultAssay(ss2) <- "ATAC_macs3"

ss2 <- RegionStats(ss2, genome = BSgenome.Drerio.UCSC.danRer11)

# link peaks to genes (nps subset)
ss2 <- LinkPeaks(
  object = ss2,
  peak.assay = "ATAC_macs3",
  expression.assay = "RNA"

)

# save Seurat objects -----------------------------------------------------

saveRDS(ss2,paste0("./data/rds/Step6_nps_var",r.variable,".rds"))

