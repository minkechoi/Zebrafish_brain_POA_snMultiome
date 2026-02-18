#!/usr/bin/Rscript

library(Seurat)
library(Signac)

set.seed(1234)

r.variable=4000
vs="06_2025_v4"

#load rds
ft_comb_seurat=readRDS(paste0("./data/rds/step3_norm_harmony_wnn_RNA_ft",r.variable,".rds"))

#fragment file path correction in case  
unix_frag=c("../10x_counts/11254_1control_pre/atac_fragments.tsv.gz",
            "../10x_counts/11254_2bPAC_pre/atac_fragments.tsv.gz",
            "../10x_counts/11254_3control_post/atac_fragments.tsv.gz",
            "../10x_counts/11254_4bPAC_post/atac_fragments.tsv.gz")

ft_comb_seurat@assays[["ATAC"]]@fragments[[1]]@path=unix_frag[1]
ft_comb_seurat@assays[["ATAC"]]@fragments[[2]]@path=unix_frag[3]
ft_comb_seurat@assays[["ATAC"]]@fragments[[3]]@path=unix_frag[2]
ft_comb_seurat@assays[["ATAC"]]@fragments[[4]]@path=unix_frag[4]

DefaultAssay(ft_comb_seurat) = "ATAC"
macs3_path = "./data/macs3_callpeak2"
dir.create(macs3_path)

macs3peaks=CallPeaks(
  object=ft_comb_seurat,
  assay="ATAC",
  group.by="wnn_cluster",
  macs2.path = "/home/mchoi/miniconda3/envs/MACS3/bin/macs3",
  outdir = macs3_path,
  fragment.tempdir = macs3_path, 
  extsize = 200,
  shift = -100,
  effective.genome.size = 1.4e+09,
  additional.args = "--bdg --SPMR -f BEDPE",
  name = sample,
  cleanup = FALSE,
  verbose = TRUE
)
saveRDS(macs3peaks, file="./data/rds/macs3_peaks2.rds")
