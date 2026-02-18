#!/usr/bin/env Rscript
#Author: Min K Choi (m.choi@exeter.ac.uk)
###packages and functions 
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
library(Seurat)
library(Signac)
library(Matrix)
library(tidyverse)
library(biomaRt)
library(org.Dr.eg.db)
library(BSgenome.Drerio.UCSC.danRer11)
library(GenomicRanges)
library(GenomicFeatures)
library(cowplot)
library(patchwork)
library(RColorBrewer)
library(viridis)
library(GenomeInfoDb)


options("scipen" = 100)
vs="06_2025_v4"
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

#functions adopted and modified from Alberto Perez-Posada @apposada
source("./ext_code/r_code/functions/sourcefolder.R")
sourceFolder("./ext_code/r_code/functions/")

# Assign arguments and paths ----------------------------------------------

scMulome_counts_DIR=paste0("../10x_counts/")
genotype=c("cont","bPAC")
condition=c("preLD","postLD")
sample_key=paste0(rep(genotype,2),"_",c("pre","pre","post","post"))
samples = c("11254_1control_pre","11254_2bPAC_pre",
            "11254_3control_post","11254_4bPAC_post")

# annotation_for Danio rerio ----------------------------------------------

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


# creating Seurat objects -------------------------------------------------
# Set multicores for matrix count
set.seed(1234)
plan("multicore", workers = 12)
options(future.globals.maxSize = 120000 * 1024^2)

for (i in 1:4) {
  #file path
  count_mtx_path = paste0(scMulome_counts_DIR,samples[i], "/filtered_feature_bc_matrix.h5")
  frag_path = paste0(scMulome_counts_DIR,samples[i], "/atac_fragments.tsv.gz")
  metadata_path = paste0(scMulome_counts_DIR,samples[i], "/per_barcode_metrics.csv")

  # Ensure required files exist
  if (!file.exists(count_mtx_path)) {
    stop("Count matrix file does not exist: ", count_mtx_path, "\n")
  }
  if (!file.exists(frag_path)) {
    stop("Fragments file does not exist: ", frag_path, "\n")
  }
  if (!file.exists(metadata_path)) {
    stop("Metadata file does not exist: ", metadata_path, "\n")
  }
  
  #read h5 file
  h5=Read10X_h5(count_mtx_path)
  
  ##SeuratObject Create
  s_obj = CreateSeuratObject(counts = h5$`Gene Expression`,
                                     assay = "RNA",
                                     project = sample_key[i])
  s_obj[['ATAC']] = CreateChromatinAssay(counts = h5$`Peaks`,
                                                 annotation = annotations,
                                                 fragments = frag_path,
                                                 sep = c(":", "-"),
                                                 genome = 'danRer11')
  
  
  #RNA_vs change
  s_obj[["RNA"]]=as(object = s_obj[["RNA"]], Class = "Assay")
  
  #normalization
  #s_obj@assays$RNA@data <- do_pf(log1p(do_pf(s_obj@assays$RNA@data)))
  
  s_obj <- FindVariableFeatures(s_obj, selection.method = "vst", nfeatures = 7000)
  s_obj[["library"]] <- sample_key[i]
  s_obj@meta.data$library <- sample_key[i]
  
  metadata= read.csv(
    file = paste(metadata_path),
    header = TRUE,
    row.names = 1
  )
  
  meta_key = paste0(sample_key[i],"_metadata")
  assign(sample_key[i],s_obj)
  assign(meta_key,metadata)
  s_obj=NULL
  print(paste0(sample_key[i],"_done."))
}

# merging sObjs -----------------------------------------------------------

add.cell.ids <- sample_key[c(1,3,2,4)]
#merge

combined_seurat = merge(get(sample_key[1]),list(get(sample_key[3]),
                                                get(sample_key[2]),
                                                get(sample_key[4])),
                        add.cell.ids = add.cell.ids, merge.data = FALSE)

Idents(combined_seurat) <- "library"

#simple check
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

# normalized merged data
#combined_seurat@assays$RNA@data <-
#  do_pf(log1p(do_pf(combined_seurat@assays$RNA@data)))

##scATAC_peak 

peaks <-  GenomicRanges::reduce(unlist(as(c(get(sample_key[1])@assays$ATAC@ranges,
                            get(sample_key[3])@assays$ATAC@ranges,
                            get(sample_key[2])@assays$ATAC@ranges,
                            get(sample_key[4])@assays$ATAC@ranges),
                          "GRangesList")))
peakwidths <- width(peaks)
peaks <- peaks[peakwidths < 10000 & peakwidths > 20]

counts_atac_merged <- FeatureMatrix(combined_seurat@assays$ATAC@fragments,
                                    features = peaks,
                                    cells = colnames(combined_seurat))
combined_seurat[['ATAC']] <- CreateChromatinAssay(counts_atac_merged,
                                                  fragments = combined_seurat@assays$ATAC@fragments,
                                                  annotation = combined_seurat@assays$ATAC@annotation,
                                                  sep = c(":","-"),
                                                  genome = "danRer11")


combined_seurat <- PercentageFeatureSet(combined_seurat, pattern = "^mt-", col.name = "percent.mt", assay = "RNA")
combined_seurat <- NucleosomeSignal(combined_seurat, assay = "ATAC")
combined_seurat <- TSSEnrichment(combined_seurat, assay = "ATAC")

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

#saveRDS(seurat_cont1,"./data/oxt_seurat_cont1.rds")
#saveRDS(seurat_cont2,"./data/oxt_seurat_cont2.rds")
#saveRDS(seurat_bPAC1,"./data/oxt_seurat_bPAC1.rds")
#saveRDS(seurat_bPAC2,"./data/oxt_seurat_bPAC2.rds")
saveRDS(combined_seurat,"./data/rds/combined_seurat_pfnorm.rds")

combined_seurat=readRDS(paste0("./data/rds/combined_seurat.rds"))

#test
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
# normalized merged data
#combined_seurat@assays$RNA@data <-
#  do_pf(log1p(do_pf(combined_seurat@assays$RNA@data)))

#cutoff_ low Quality cells.

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


saveRDS(ft_comb_seurat,paste0("./data/rds/ft_comb_seurat.rds"))
combined_seurat=NULL
gc()


