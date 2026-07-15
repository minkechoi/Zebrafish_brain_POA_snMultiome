#!/usr/bin/env Rscript
# =============================================================================
# 01_seurat_objects_wo_pf.R
# -----------------------------------------------------------------------------
# Purpose : Build the primary single-nucleus multiome (snRNA + snATAC) Seurat
#           object for the zebrafish preoptic area (POA). Reads Cell Ranger ARC
#           per-sample outputs, creates one Seurat object per sample (paired
#           Gene Expression + ATAC peaks), merges the four libraries, recomputes
#           a common ATAC peak set, computes QC metrics, and filters out
#           low-quality nuclei.
# Design  : 2 genotypes (control, bPAC) x 2 conditions (preLD, postLD) = 4 libs.
#           "pf" in the filename refers to proportional-fitting normalisation,
#           which is intentionally NOT applied here (raw counts are retained;
#           see commented-out do_pf() calls).
# Inputs  : ../10x_counts/<sample>/filtered_feature_bc_matrix.h5
#           ../10x_counts/<sample>/atac_fragments.tsv.gz
#           ../10x_counts/<sample>/per_barcode_metrics.csv
#           AnnotationHub EnsDb (Danio rerio, Ensembl 111) for gene models.
# Outputs : ./<vs>/data/rds/combined_seurat_pfnorm.rds  (merged, pre-filter)
#           ./<vs>/data/rds/ft_comb_seurat.rds          (QC-filtered)
#           QC violin/scatter/histogram plots under ./<vs>/figures/
# Author  : Min K Choi (m.choi@exeter.ac.uk)
# =============================================================================

### packages and functions
rm(list = ls(all.names = TRUE)) # clear all objects, including hidden ones, for a clean session
gc()                            # free memory and report current memory usage
library(Seurat)                 # core single-cell object model and workflow
library(Signac)                 # scATAC-seq extension for Seurat (ChromatinAssay)
library(Matrix)                 # sparse matrix support for count matrices
library(tidyverse)              # data wrangling / plotting (dplyr, ggplot2, ...)
library(biomaRt)                # Ensembl BioMart gene annotation queries
library(org.Dr.eg.db)           # zebrafish (Danio rerio) gene ID mappings
library(BSgenome.Drerio.UCSC.danRer11) # danRer11 genome sequence
library(GenomicRanges)          # genomic interval arithmetic
library(GenomicFeatures)        # transcript/annotation database utilities
library(cowplot)                # figure composition
library(patchwork)              # combine ggplots with | and /
library(RColorBrewer)           # colour palettes
library(viridis)                # perceptually uniform colour palettes
library(GenomeInfoDb)           # seqinfo / chromosome-style helpers


options("scipen" = 100)         # prefer fixed over scientific notation in output
vs="06_2025_v4"                 # version tag; all outputs are written under ./<vs>/
# Create the standard project sub-folder layout for this version
dir.create(paste0("./",vs))
dir.create(paste0("./",vs,"/data"))
dir.create(paste0("./",vs,"/figures"))
dir.create(paste0("./",vs,"/outputs"))
dir.create(paste0("./",vs,"/codes"))


#set working dir
setwd(paste0("./",vs))


#save session info and Rstudio vs info for reproducibility
writeLines(capture.output(sessionInfo()), "./sessionInfo.txt")
#writeLines(capture.output(rstudioapi::vsInfo()), "./vsInfo.txt")

# Load helper functions (sourceFolder() sources every .R in a directory).
# functions adopted and modified from Alberto Perez-Posada @apposada
source("./ext_code/r_code/functions/sourcefolder.R")
sourceFolder("./ext_code/r_code/functions/")

# Assign arguments and paths ----------------------------------------------
# Define the four libraries and the mapping to their Cell Ranger ARC output dirs.

scMulome_counts_DIR=paste0("../10x_counts/")     # root of Cell Ranger ARC counts
genotype=c("cont","bPAC")                        # two genotypes
condition=c("preLD","postLD")                    # two light/dark conditions
sample_key=paste0(rep(genotype,2),"_",c("pre","pre","post","post")) # short sample labels
samples = c("11254_1control_pre","11254_2bPAC_pre",
            "11254_3control_post","11254_4bPAC_post")               # Cell Ranger folder names

# annotation_for Danio rerio ----------------------------------------------
# Fetch the Ensembl 111 gene models for zebrafish and convert to UCSC/danRer11
# style so they match the ATAC fragment coordinates.

library(AnnotationHub)
#library(biomaRt)
#library(org.Dr.eg.db)

ah = AnnotationHub()
ensdbs <- query(ah, c("Danio rerio"))                                  # all D. rerio EnsDb records
ensdb_id <- ensdbs$ah_id[grep(paste0(" 111 EnsDb"), ensdbs$title)]     # pick Ensembl release 111
ensdb <- ensdbs[[ensdb_id]]
annotations <- GetGRangesFromEnsDb(ensdb = ensdb, standard.chromosomes =T) # gene models as GRanges
seqlevelsStyle(annotations) = "UCSC"                                   # chr-style to match ATAC
genome(annotations) = "danRer11"


# creating Seurat objects -------------------------------------------------
# For each of the 4 samples: read the multiome H5, build a Seurat object with a
# paired RNA assay and ATAC ChromatinAssay, and attach per-barcode metadata.
# Set multicores for matrix count
set.seed(1234)                                       # reproducible randomness
plan("multicore", workers = 12)                      # parallel backend for heavy steps
options(future.globals.maxSize = 120000 * 1024^2)    # raise memory cap for large objects

for (i in 1:4) {
  #file path
  count_mtx_path = paste0(scMulome_counts_DIR,samples[i], "/filtered_feature_bc_matrix.h5")
  frag_path = paste0(scMulome_counts_DIR,samples[i], "/atac_fragments.tsv.gz")
  metadata_path = paste0(scMulome_counts_DIR,samples[i], "/per_barcode_metrics.csv")

  # Ensure required files exist (fail fast with an informative message)
  if (!file.exists(count_mtx_path)) {
    stop("Count matrix file does not exist: ", count_mtx_path, "\n")
  }
  if (!file.exists(frag_path)) {
    stop("Fragments file does not exist: ", frag_path, "\n")
  }
  if (!file.exists(metadata_path)) {
    stop("Metadata file does not exist: ", metadata_path, "\n")
  }

  #read h5 file (contains both "Gene Expression" and "Peaks" matrices)
  h5=Read10X_h5(count_mtx_path)

  ##SeuratObject Create (RNA assay from the Gene Expression matrix)
  s_obj = CreateSeuratObject(counts = h5$`Gene Expression`,
                                     assay = "RNA",
                                     project = sample_key[i])
  # Add the ATAC peak matrix as a Signac ChromatinAssay (with gene annotation + fragments)
  s_obj[['ATAC']] = CreateChromatinAssay(counts = h5$`Peaks`,
                                                 annotation = annotations,
                                                 fragments = frag_path,
                                                 sep = c(":", "-"),
                                                 genome = 'danRer11')


  #RNA_vs change: coerce RNA to the classic "Assay" class (v3-style) for downstream compatibility
  s_obj[["RNA"]]=as(object = s_obj[["RNA"]], Class = "Assay")

  #normalization (proportional-fitting normalisation intentionally disabled here)
  #s_obj@assays$RNA@data <- do_pf(log1p(do_pf(s_obj@assays$RNA@data)))

  s_obj <- FindVariableFeatures(s_obj, selection.method = "vst", nfeatures = 7000) # top 7000 HVGs
  s_obj[["library"]] <- sample_key[i]              # tag each cell with its library
  s_obj@meta.data$library <- sample_key[i]

  metadata= read.csv(
    file = paste(metadata_path),
    header = TRUE,
    row.names = 1
  )

  # Store the object and its metadata under names derived from sample_key
  meta_key = paste0(sample_key[i],"_metadata")
  assign(sample_key[i],s_obj)
  assign(meta_key,metadata)
  s_obj=NULL                                        # release memory before next iteration
  print(paste0(sample_key[i],"_done."))
}

# merging sObjs -----------------------------------------------------------
# Merge the four per-sample objects into one. Cell barcodes are prefixed with
# the library id to keep them unique across samples.

add.cell.ids <- sample_key[c(1,3,2,4)]              # prefix order matches merge order below
#merge

combined_seurat = merge(get(sample_key[1]),list(get(sample_key[3]),
                                                get(sample_key[2]),
                                                get(sample_key[4])),
                        add.cell.ids = add.cell.ids, merge.data = FALSE)

Idents(combined_seurat) <- "library"               # set library as the active identity

#simple check: compare gene detection across libraries
DefaultAssay(combined_seurat)="RNA"
ComparisonLibrary_nfeat <- VlnPlot(combined_seurat,
                                   assay = "RNA",
                                   features = "nFeature_RNA",
                                   group.by = "orig.ident")+
  ggtitle("_combined_seurat_merge_nFeature_By_Library_highgene_cutoff.pdf")
print(ComparisonLibrary_nfeat)

ComparisonLibrary_ncount <- VlnPlot(combined_seurat, "nCount_RNA", group.by = "orig.ident")+
  ggtitle("_combined_seurat_merge_nCount_By_Library_highgene_cutoff.pdf")
print(ComparisonLibrary_ncount)

# normalized merged data (proportional-fitting normalisation intentionally disabled)
#combined_seurat@assays$RNA@data <-
#  do_pf(log1p(do_pf(combined_seurat@assays$RNA@data)))

##scATAC_peak
# Per-sample ATAC peaks were called independently, so build one unified peak set:
# take the union of all sample peak ranges and reduce (merge) overlaps.

peaks <-  GenomicRanges::reduce(unlist(as(c(get(sample_key[1])@assays$ATAC@ranges,
                            get(sample_key[3])@assays$ATAC@ranges,
                            get(sample_key[2])@assays$ATAC@ranges,
                            get(sample_key[4])@assays$ATAC@ranges),
                          "GRangesList")))
peakwidths <- width(peaks)
peaks <- peaks[peakwidths < 10000 & peakwidths > 20]  # drop implausibly wide/narrow peaks

# Re-quantify the merged peak set across all cells and rebuild the ATAC assay
counts_atac_merged <- FeatureMatrix(combined_seurat@assays$ATAC@fragments,
                                    features = peaks,
                                    cells = colnames(combined_seurat))
combined_seurat[['ATAC']] <- CreateChromatinAssay(counts_atac_merged,
                                                  fragments = combined_seurat@assays$ATAC@fragments,
                                                  annotation = combined_seurat@assays$ATAC@annotation,
                                                  sep = c(":","-"),
                                                  genome = "danRer11")


# QC metrics: mitochondrial %, nucleosome banding, and TSS enrichment
combined_seurat <- PercentageFeatureSet(combined_seurat, pattern = "^mt-", col.name = "percent.mt", assay = "RNA")
combined_seurat <- NucleosomeSignal(combined_seurat, assay = "ATAC")
combined_seurat <- TSSEnrichment(combined_seurat, assay = "ATAC")

# QC overview plots (pre-filter)
a=VlnPlot(combined_seurat,
          features = c("nFeature_RNA",
                       "percent.mt",
                       "nFeature_ATAC",
                       "TSS.enrichment",
                       "nucleosome_signal"),
          ncol = 5,
          pt.size = 0)

b <- FeatureScatter(combined_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


tiff("./figures/qc_violin_plots_prefilter.tiff",width = 40,height = 15,units = "cm",
     res = 300,compression = "lzw")
print(a|b)+ patchwork::plot_layout(widths = c(4,1))
dev.off()

dir.create(paste0("./data/rds"))

# (per-sample saves left disabled; only the merged object is persisted)
#saveRDS(seurat_cont1,"./data/oxt_seurat_cont1.rds")
#saveRDS(seurat_cont2,"./data/oxt_seurat_cont2.rds")
#saveRDS(seurat_bPAC1,"./data/oxt_seurat_bPAC1.rds")
#saveRDS(seurat_bPAC2,"./data/oxt_seurat_bPAC2.rds")
saveRDS(combined_seurat,"./data/rds/combined_seurat_pfnorm.rds")

combined_seurat=readRDS(paste0("./data/rds/combined_seurat.rds"))

# test: distributions of QC metrics on a log10 scale to help choose cutoffs
pdf(paste0("./figures/nCounts_nFeature_distribution_histPlot.pdf", width=9,height=6))
par(mfrow=c(2,4))
hist(log10(combined_seurat$nCount_RNA),n=1000,col="darkgrey",border="darkgrey",main="nCount_RNA")
hist(log10(combined_seurat$nFeature_RNA),n=1000,col="darkgrey",border="darkgrey",main="nFeature_RNA")
hist(log10(combined_seurat$nCount_ATAC),n=1000,col="darkgrey",border="darkgrey",main="nCount_ATAC")
hist(log10(combined_seurat$nFeature_ATAC),n=1000,col="darkgrey",border="darkgrey",main="nFeature_ATAC")
hist(combined_seurat$nucleosome_signal,n=1000,col="darkgrey",border="darkgrey",main="Nucleosome_signal")
hist(combined_seurat$TSS.enrichment,n=1000,col="darkgrey",border="darkgrey",main="TSS_enrichment")
hist(combined_seurat$percent.mt,n=1000,col="darkgrey",border="darkgrey",main="Pct_mt")
dev.off()


##filtering, sub-setting

# subset
# normalized merged data (proportional-fitting normalisation intentionally disabled)
#combined_seurat@assays$RNA@data <-
#  do_pf(log1p(do_pf(combined_seurat@assays$RNA@data)))

# cutoff_ low Quality cells.
# Keep nuclei within sensible RNA depth/complexity bounds, low mito contamination,
# and good ATAC quality (TSS enrichment high, nucleosome signal low).

ft_comb_seurat= subset(x = combined_seurat,
                       nCount_RNA >= 200 &
                         nCount_RNA <= 25000 &
                         nFeature_RNA >= 400&
                         nFeature_RNA <= 7000 &
                         percent.mt < 2 &
                         #nFeature_ATAC > as.numeric(quantile(ss$nCount_ATAC, probs = 0.05)) &
                         TSS.enrichment > 2 &
                         nucleosome_signal < 2
)

gc()

# QC overview plots (post-filter) to confirm the cutoffs cleaned the data
a=VlnPlot(ft_comb_seurat,
          features = c("nFeature_RNA",
                       "percent.mt",
                       "nFeature_ATAC",
                       "TSS.enrichment",
                       "nucleosome_signal"),
          ncol = 5,
          pt.size = 0)

b <- FeatureScatter(ft_comb_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

tiff("./figures/qc_violin_plots_postfilter.tiff",width = 40,height = 15,units = "cm",
     res = 300,compression = "lzw")
print(a|b)+ patchwork::plot_layout(widths = c(4,1))
dev.off()


saveRDS(ft_comb_seurat,paste0("./data/rds/ft_comb_seurat.rds")) # QC-filtered object for step 03
combined_seurat=NULL
gc()


