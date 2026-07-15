#!/usr/bin/Rscript
# =============================================================================
# 03_2_macs3_unix_run.R
# -----------------------------------------------------------------------------
# Purpose : Call ATAC peaks per cluster with MACS3, run on a Unix/HPC machine.
#           Loads the normalised WNN object, repoints ATAC fragment paths to the
#           Unix filesystem, then calls cluster-resolved peaks via Signac's
#           CallPeaks() wrapper around the MACS3 binary.
# Inputs  : ./data/rds/step3_norm_harmony_wnn_RNA_ft<r.variable>.rds
#           ../10x_counts/<sample>/atac_fragments.tsv.gz
# Output  : ./data/rds/macs3_peaks2.rds  (+ raw MACS3 output in macs3_callpeak2/)
# Note    : Grouped by `wnn_cluster`; uses ATAC-appropriate MACS3 settings
#           (BEDPE input, shift/extsize for Tn5, zebrafish effective genome).
# =============================================================================

library(Seurat)
library(Signac)

set.seed(1234)

r.variable=4000       # number of variable features used upstream (encoded in the RDS name)
vs="06_2025_v4"       # version tag

#load rds (normalised RNA+ATAC WNN object from step 03_0)
ft_comb_seurat=readRDS(paste0("./data/rds/step3_norm_harmony_wnn_RNA_ft",r.variable,".rds"))

# fragment file path correction in case the object was created on a different (e.g. Windows) machine
unix_frag=c("../10x_counts/11254_1control_pre/atac_fragments.tsv.gz",
            "../10x_counts/11254_2bPAC_pre/atac_fragments.tsv.gz",
            "../10x_counts/11254_3control_post/atac_fragments.tsv.gz",
            "../10x_counts/11254_4bPAC_post/atac_fragments.tsv.gz")

# Reassign each fragment object's path (note: index order maps stored order -> sample order)
ft_comb_seurat@assays[["ATAC"]]@fragments[[1]]@path=unix_frag[1]
ft_comb_seurat@assays[["ATAC"]]@fragments[[2]]@path=unix_frag[3]
ft_comb_seurat@assays[["ATAC"]]@fragments[[3]]@path=unix_frag[2]
ft_comb_seurat@assays[["ATAC"]]@fragments[[4]]@path=unix_frag[4]

DefaultAssay(ft_comb_seurat) = "ATAC"
macs3_path = "./data/macs3_callpeak2"
dir.create(macs3_path)

# Per-cluster peak calling with MACS3 (Tn5-shifted, BEDPE, single-per-million-reads)
macs3peaks=CallPeaks(
  object=ft_comb_seurat,
  assay="ATAC",
  group.by="wnn_cluster",               # call peaks separately within each WNN cluster
  macs2.path = "/home/mchoi/miniconda3/envs/MACS3/bin/macs3",
  outdir = macs3_path,
  fragment.tempdir = macs3_path,
  extsize = 200,                          # fragment extension size
  shift = -100,                           # shift reads to centre on the Tn5 cut site
  effective.genome.size = 1.4e+09,        # danRer11 mappable genome size
  additional.args = "--bdg --SPMR -f BEDPE", # bedGraph output, signal-per-million-reads, paired-end
  name = sample,
  cleanup = FALSE,
  verbose = TRUE
)
saveRDS(macs3peaks, file="./data/rds/macs3_peaks2.rds")
