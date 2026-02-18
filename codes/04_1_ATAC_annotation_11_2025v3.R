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

#set working dir
if(getwd()!=paste0("D:/projects/scMultiome_oxt/",vs)){
  setwd(paste0("./",vs))
}

#functions adopted and modified from Alberto Perez-Posada @apposada
source("./ext_code/r_code/functions/sourcefolder.R")
sourceFolder("./ext_code/r_code/functions/")

#other temp function
'%!in%' <- function(x,y)!('%in%'(x,y))
#load data
#ss1=readRDS(paste0("./data/rds/step4_var",r.variable,".rds"))

# ATAC_link ---------------------------------------------------------------

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

ss1 <- RegionStats(ss1, genome = BSgenome.Drerio.UCSC.danRer11)

# link peaks to genes
ss1 <- LinkPeaks(
  object = ss1,
  peak.assay = "ATAC_macs3",
  expression.assay = "RNA",
    
)

# save Seurat objects -----------------------------------------------------

saveRDS(ss1,paste0("./data/rds/Step4_var",r.variable,".rds"))




