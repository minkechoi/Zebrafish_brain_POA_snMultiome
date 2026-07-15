#!/usr/bin/env Rscript
# =============================================================================
# 03_1_macs3_peak.R
# -----------------------------------------------------------------------------
# Purpose : Replace the default (Cell Ranger) ATAC peak set with a cluster-aware
#           MACS3 peak set. Drives the MACS3 run (03_2), filters peaks (standard
#           chromosomes + blacklist), re-quantifies them into a new ChromatinAssay,
#           computes FRiP QC, and swaps it in as the working ATAC assay.
# Inputs  : ./data/rds/step3_norm_harmony_wnn_RNA_ft<r.variable>.rds
#           ./data/rds/macs3_peaks2.rds (from 03_2)
#           ./data/Blacklist_danRer10_to_danRer11_YueLab_srt.bed
# Output  : ./data/rds/step3_norm_harmony_wnn_RNA_ft_macs<r.variable>.rds
# Note    : MACS3 must be run on Unix/Mac (not Windows); see 03_2 for the call.
# =============================================================================

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc()

#######

library(Seurat)
library(Signac)
library(tidyverse)
library(BSgenome.Drerio.UCSC.danRer11)   # danRer11 genome sequence
library(GenomicRanges)
library(GenomicFeatures)
library(GenomeInfoDb)
library(rtracklayer)                     # import BED (blacklist) files
library(future)                          # parallel backend
library(viridis)

set.seed(1234)
#plan("multicore", workers = 8)
#options(future.globals.maxSize = 120000 * 1024^2)
options(future.globals.maxSize = 8000 * 1024^2)   # memory cap for future globals

r.variable=4000
vs="06_2025_v4"
setwd(paste0("./",vs))

# MACS_peak calling -------------------------------------------------------

###MACS is not available in windows environment.Use unix or Mac environment
####### 1. Macs peak calling using filtered cells #######
# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth

###### Run MACS in the unix or mac env which MACS3 pre-installed  !!!!!!!

# Run the actual MACS3 peak calling (03_2) — produces ./data/rds/macs3_peaks2.rds
source("./codes/03_2_macs3_unix_run.R")
#or
#system("bash ./codes/sh_files/03_2_macs3_unix_run.sh")
#continue back to windows environment if your are working in windows env.

#######

##load s.obj and macs results
ft_comb_seurat=readRDS(paste0("./data/rds/step3_norm_harmony_wnn_RNA_ft",r.variable,".rds"))
DefaultAssay(ft_comb_seurat)="ATAC"
macs3peaks=readRDS(file="./data/rds/macs3_peaks2.rds")

length(macs3peaks)                        # number of raw MACS3 peaks
write.table(as.data.frame(macs3peaks),file="./outputs/macs3_peaks2.txt",sep="\t",quote=FALSE,row.names=FALSE)

# remove peaks on nonstandard chromosomes (we dont have, since we removed the nonstandard chr during cellranger-arc count)
macs3peaks = keepStandardChromosomes(macs3peaks, pruning.mode = "coarse")

# remove peaks in genomic blacklist regions (ENCODE-style artefact regions)
blacklist_danRer11 =rtracklayer::import("./data/Blacklist_danRer10_to_danRer11_YueLab_srt.bed")
macs3peaks = subsetByOverlaps(x = macs3peaks, ranges = blacklist_danRer11, invert = TRUE)

saveRDS(macs3peaks, file="./data/rds/macs3_peaks_afterfiltering2.rds")
#macs3peaks= readRDS("./data/rds/macs3_peaks_afterfiltering.rds")

# Gene annotation (UCSC/danRer11) for the new peak assay
library(AnnotationHub)

ah = AnnotationHub()
ensdbs <- query(ah, c("Danio rerio"))
ensdb_id <- ensdbs$ah_id[grep(paste0(" 111 EnsDb"), ensdbs$title)]
ensdb <- ensdbs[[ensdb_id]]
annotations <- GetGRangesFromEnsDb(ensdb = ensdb, standard.chromosomes =T)
seqlevelsStyle(annotations) = "UCSC"
genome(annotations) = "danRer11"


# create macs3peaks-cell matrix (count fragments per MACS3 peak per cell)
ft_comb_seurat_macs3counts = FeatureMatrix(
  fragments = Fragments(ft_comb_seurat),
  features = macs3peaks,
  cells = colnames(ft_comb_seurat),
  process_n = 2000
)

# create a new assay using the macs3 peak set and add it to the Seurat object
ft_comb_seurat[["ATAC_cluster"]] = CreateChromatinAssay(
  counts = ft_comb_seurat_macs3counts,
  fragments = Fragments(ft_comb_seurat),
  annotation = annotations,
  min.cells = 5
)

# Calculate reads in peak ratio (FRiP): a per-cell ATAC quality metric
DefaultAssay(ft_comb_seurat)='ATAC_cluster'
ft_comb_seurat = FRiP(object = ft_comb_seurat,
                      assay = 'ATAC_cluster',
                      total.fragments = 'nCount_ATAC_cluster')

#FRIP (compute manually as well and store as FRiP2)
peak.data <- GetAssayData(object = ft_comb_seurat, assay = "ATAC_cluster",
                          layer = "counts")
total_fragments_cell <- ft_comb_seurat[[]][["nCount_ATAC_cluster"]]
peak.counts <- colSums(x = peak.data)
frip <- peak.counts/total_fragments_cell
ft_comb_seurat <- AddMetaData(object = ft_comb_seurat, metadata = frip, col.name = "FRiP2")


# add number of peak region fragments for Signac called peaks (equals the value in the "nCount_peaks" value in the object)
ft_comb_seurat$macs3peaks_region_fragments2 =
  ft_comb_seurat$FRiP2 * ft_comb_seurat$nCount_ATAC_cluster

#  Add pct_reads_in_peaks for Signac called peaks to seurat object
ft_comb_seurat$pct_reads_in_macs3peaks2 = ft_comb_seurat$macs3peaks_region_fragments2 / ft_comb_seurat$nCount_ATAC_cluster * 100

#replace ATAC to ATAC_cluster (make the MACS3 peak set the primary ATAC assay)

ft_comb_seurat[["ATAC"]] = ft_comb_seurat[["ATAC_cluster"]]

DefaultAssay(ft_comb_seurat)='ATAC'
ft_comb_seurat[["ATAC_cluster"]]=NULL           # drop the temporary duplicate assay
saveRDS(ft_comb_seurat,paste0("./data/rds/step3_norm_harmony_wnn_RNA_ft_macs",r.variable,".rds"))

