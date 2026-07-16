https://doi.org/10.5281/zenodo.21386287

# Zebrafish hypothalamic preoptic area (POA) analysis

Single-nuclei Multiome-sequencing data analysis (snRNA-seq + snATAC-seq).

Code and data for the manuscript *"Latent priming of hypothalamic neuroendocrine populations links developmental glucocorticoid excess to adult stress hypersensitivity"*.

---

## About

This repository hosts the code used to perform the analyses from the manuscript:

> Min-Kyeung Choi, Anna Tochwin, Alberto Perez Posada, Jordi Solana, Soojin Ryu. *Latent priming of hypothalamic neuroendocrine populations links developmental glucocorticoid excess to adult stress hypersensitivity*, 2026.

Here you can find all the code that was used to generate the panels for the majority of the main and supplementary figures, as well as the majority of the supplementary files.

---

## Data availability

**Sequencing data and processed files**
GEO: **GSE320318** — snMultiome studies in zebrafish hypothalamic preoptic area following developmental glucocorticoid overexposure. Currently available for reviewers only.

**Confocal microscopy images**
Original confocal microscopy image files will be available on figshare (`10.6084/m9.figshare.32719890`). Currently available for reviewers only.

---

## Background

Excess stress exposure during early life can affect sensitivity to subsequent stress both in humans and animal models. Studies show that glucocorticoids (GCs) play a pivotal role in mediating the long-lasting effects of stress. The hypothalamus represents an important regulatory brain region for stress hormone control. However, it is currently not known whether and how developmental GC over-exposure (dGC-OE) alters hypothalamic cells to shape their adulthood response to stress. To tackle this knowledge gap, we developed a double-hit stress zebrafish model, which combines developmental GC over-exposure (dGC-OE) and acute predatory-like stress exposure by looming dots (LD) in adulthood. The preoptic area of zebrafish — the zebrafish equivalent of the mammalian anterior hypothalamic area including the paraventricular hypothalamic nucleus (PVN) — was analyzed using snMultiome-seq (GEX + ATAC).

---

## 10X Multiome-seq

Single nuclei were isolated from snap-frozen preoptic area tissues of 8-month-old negative control and transgenic (dGC-OE group, `Tg(star:bPAC-2A-tdTomato)ue300+/-`) zebrafish, with and without adulthood acute stress exposure. Tissue samples were collected 15 minutes post-acute stress and pooled (n = 12 per sample) to ensure sufficient material for analysis. Joint snRNA- and snATAC libraries were prepared using the 10x Genomics Chromium Next GEM Single Cell Multiome ATAC + Gene Expression platform.

The experimental design is a 2 × 2 layout — **2 genotypes** (control, bPAC / dGC-OE) × **2 conditions** (pre-LD, post-LD) = 4 libraries.

---

## Basic structure of the repository

```
Zebrafish_brain_POA_snMultiome/
├── README.md                 # This file
├── LICENSE
├── sessionInfo.txt           # R session / package versions for reproducibility
├── codes/                    # All analysis scripts (R), numbered by pipeline stage
│   ├── Figure/               # Scripts for specific manuscript figure panels (+ Cytoscape session)
│   └── sh_files/             # Shell scripts (Cell Ranger ARC, MACS3, HOMER) run on Unix/HPC
├── data/                     # Small reference files used by the scripts
│   ├── Blacklist_danRer10_to_danRer11_YueLab_srt.bed   # ATAC blacklist regions (danRer11)
│   └── danRer11.gimme.vertebrate.v5.0.motif2factors.txt # motif → TF mapping (ANANSE/GimmeMotifs)
└── ext_code/                 # External helper code (see credits below)
    ├── comparABle-main/      # comparABle utilities (tidyup, ensemble clustering)
    └── r_code/               # shared R helper functions (sourced by the scripts)
```

Large inputs (Cell Ranger ARC outputs, per-sample `.h5` / `atac_fragments.tsv.gz`) and generated outputs (`.rds` objects, figures, tables) are **not** stored in the repository; they are produced by the scripts and/or available from GEO (see *Data availability*).

### Analysis scripts (`codes/`)

The R scripts are numbered to reflect the order of the pipeline. Each script begins with a header describing its purpose, inputs, and outputs.

| Stage | Script(s) | Description |
|-------|-----------|-------------|
| **01 — Object build & QC** | `01_seurat_objects_wo_pf.R` | Read Cell Ranger ARC outputs, build the paired RNA + ATAC Seurat object, merge the 4 libraries, compute QC metrics, filter low-quality nuclei. |
| **03 — Normalisation & peaks** | `03_0_seurat_objects_normalization2_*.R` | SCTransform + cell-cycle regression, PCA, Harmony batch correction, WNN joint RNA+ATAC clustering/UMAP. |
| | `03_1_macs3_peak.R`, `03_2_macs3_unix_run.R` | Cluster-aware MACS3 peak calling and rebuild of the ATAC assay (FRiP QC). |
| **04 — Annotation & scoring** | `04_0_cluster_annotationv4_*.R` | Cluster annotation into named cell types (reference enrichment + "annotation voting"). |
| | `04_1_ATAC_annotation_*.R`, `06_ATAC_annotation_*.R` | Peak-to-gene linkage (Signac `LinkPeaks`). |
| | `04_a_AUC_UCell_score_calculation.R` | Stress-signature scoring (IEGs / GRs / fkbp5) with AUCell and UCell. |
| | `04_b_MiloR_run.R` | Differential abundance (cell-composition) testing with miloR. |
| | `04_c_hclust_cluster_06_2025.R`, `04_c_hclust_cluster_06_2026.R` | Hierarchical cell-type tree from pseudobulk co-occurrence clustering. |
| | `scType.R` | Auxiliary marker-based cell-type annotation (ScType). |
| **05 — Cluster selection** | `05_cluster_selection_*.R` | Select light/dark (LD) stress-responsive cell types from signature scores. |
| | `05_a_nps_sub_umaps.R` | Focused analysis of neuropeptidergic (neuroendocrine) neurons. |
| **07 — DEG / DAR** | `07_DEG_DAR_pseudobulk_v2.R` | Pseudobulk differential expression (DEG) and accessibility (DAR) per cell type (DESeq2). |
| **08 — WGCNA** | `08_a_wgcna_prep_psudobilk.R` | Pseudobulk preparation for WGCNA. |
| | `08_WGCNA_v3_*.R` | WGCNA driver: modules, module–trait relations, TF connectivity, GO / motif enrichment. |
| | `08_b_WGCNA_functions.R`, `08_c_WGCNA_graph_functions.R` | WGCNA and network-graph function libraries. |
| | `08_d_DEG_modules_net_*.R` | Integrate DEGs with WGCNA modules; export Cytoscape network. |
| **09 — Gene-regulatory networks** | `09_annasi_2025.R` | Export pseudobulk RNA/ATAC for scANANSE / ANANSE GRN inference. |
| | `09_c_motif.R`, `motif.R` | TF-motif enrichment on differentially accessible regions (identical scripts). |
| | `ananse_graph_v3.R`, `ananse_graph for average_v2.R`, `ananse_graph_function.R` | Build and analyse ANANSE gene-regulatory-network graphs. |
| **Figures** | `Figure/FIG5D.R` | Figure 5 panels (TF→target heatmaps; acute-stress pathway graphs). |
| | `Figure/Fig3D_network_cytoscape.cys` | Cytoscape session for the Figure 3D network. |
| **Shell scripts** | `sh_files/*.sh` | Cell Ranger ARC reference building & counting, MACS3, and HOMER motif finding (run on Unix/HPC). |

---

## ext_code

Adapted code from **Alberto Perez-Posada** (@apposada): https://github.com/scbe-lab/regulatory_logic
