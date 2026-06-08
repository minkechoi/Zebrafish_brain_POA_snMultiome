## Genome Visualisation
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
source("~/projects/smed_cisreg/code/r_code/functions/genomeviz.R")
options(ucscChromosomeNames=FALSE)

# Load the genomic coordinates
coords_of_interest <- read.delim2(
  file = "~/projects/smed_cisreg/outputs/genomeviz/genomic_coordinates_of_interest.tsv",
  header = TRUE
)

# Create the chromosome object
bands <- as.data.frame(readRDS("~/projects/smed_cisreg/outputs/genomeviz/cytobands.rds"))
itrack <- IdeogramTrack(genome="Smed", name = "",bands=bands)

# Add coordinate ruler
axisTrack <- GenomeAxisTrack()

# Add genemodels
smed_genemodels <- readRDS("~/projects/smed_cisreg/outputs/genomeviz/smed_cisreg_genomeviz_genemodels.rds")
smed_genemodels_longestIso <- smed_genemodels[smed_genemodels$transcript %in% rosetta$longest_isoform]
gm <- GeneRegionTrack(range=smed_genemodels_longestIso, fill = "darkblue", col = NULL)

# Add peaks scATAC
smed_peaks_scatac_bed <- import.bed("/mnt/sda/alberto/projects/smed_cisreg/outputs/scatac/NX73/outs/filtered_peak_bc_matrix/peaks.bed")
peaks_sc <- AnnotationTrack(range=smed_peaks_scatac_bed, genome="Smed", name = "peaks (scATAC)", fill = "#0c9e7e", col = NULL, background.title = "#0c9e7e")

# Add peaks bulk
smed_peaks_bed <- import.bed("/mnt/sda/alberto/colabos/virginia_atac/samples_virginia/peakcalling/peaks_virginia.bed")
peaks <- AnnotationTrack(range=smed_peaks_bed, genome="Smed", name = "peaks", fill = "#468499", col = NULL, background.title = "#468499")

# Load the track information table
tracks_table <- read.delim2(
  file = "~/projects/smed_cisreg/outputs/genomeviz/smed_cisreg_tracks.tsv",
  header = TRUE
)

tracks_table <-
  tracks_table[-c(1,2,13),]

tracks_table$ctype <- 
  factor(
    tracks_table$ctype,
    levels = c(
      "neoblast",
      "early_epidermal_progenitors",
      "late_epidermal_progenitors",
      "epidermis",
      "phagocytes",
      "basal_goblet",
      "muscle",
      "neurons",
      "parenchyma",
      "protonephridia",
      "secretory"
    )
  )

tracks_table <- tracks_table[order(tracks_table$ctype,tracks_table$data),]
tracks_table_scatac <- tracks_table[tracks_table$data=="ATAC",]
tracks_table_scrna <- tracks_table[tracks_table$data=="RNA",]

# Add smed bed file as a way to retrieve gene coordinates

smed_bed <- read.delim2("outputs/associate_peaks_genes/smed.bed", header = FALSE, col.names=c("chrom","from","to","name","dot","strand"))

## Broad RNA markers to detect OCRs nearby also open

Idents(smed_cisreg_scrna) <- smed_cisreg_scrna$broadtype

markers_broad <- 
  FindAllMarkers(
    smed_cisreg_scrna,
    only.pos = TRUE, 
    return.thresh = 1,
    logfc.threshold = 0
  )

Idents(smed_cisreg_scrna) <- smed_cisreg_scrna$ctype

markers_broad$cluster <- factor(markers_broad$cluster, levels = c(smed_ctypes$broadtype))

markers_broad_top <- markers_broad %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC)

markers_broad_coords <-
  smed_bed[smed_bed$name %in% markers_broad_top$gene,c(1:4)]

pdf("graphics/markers_broad_gviz.pdf", height = 8, width = 2)
for(i in 1:nrow(markers_broad_coords)){
  
  region_coords <- markers_broad_coords[i,]
  l <- abs(region_coords[1,2]-region_coords[1,3])
  region_coords[1,2] <- region_coords[1,2] - (l*0.5)
  region_coords[1,3] <- region_coords[1,3] + (l*0.5)
  
  # Plotting
  genomeviz(
    coords = region_coords,
    tracks = tracks_table,
    basic_tracks = list(gm, peaks, peaks_sc),
    y_axis = 100
  )
  
  message(i)
}
dev.off()

# SPECIFIC CELLTYPE CO-MARKERS AFTER MANUAL CURATION

featplot_gviz <- function(x){
  FeaturePlot(
    smed_cisreg_scrna,
    features = x,
    cols = c("#DAE7F2","#8000A8"),
    order = TRUE,
    pt.size = 1
  )+NoAxes()+NoLegend()
  }

# First the scrna
## Early Epid Prog
png("graphics/featplot_scrna_early_epid_prog.png", height = 300, width = 300)
featplot_gviz("h1SMcG0014249")
dev.off()

## Late Epid Prog
png("graphics/featplot_scrna_late_epid_prog.png", height = 300, width = 300)
featplot_gviz("h1SMnG0022130")
dev.off()

## Epid
png("graphics/featplot_scrna_epid.png", height = 300, width = 300)
featplot_gviz("h1SMcG0016267")
dev.off()

## Phag
png("graphics/featplot_scrna_phag.png", height = 300, width = 300)
featplot_gviz("h1SMcG0009317")
dev.off()

## Basal/Goblet
png("graphics/featplot_scrna_basal_goblet.png", height = 300, width = 300)
featplot_gviz("h1SMcG0002232")
dev.off()

## Musc
png("graphics/featplot_scrna_musc.png", height = 300, width = 300)
featplot_gviz("h1SMcG0003060")
dev.off()

## Neur
png("graphics/featplot_scrna_neur.png", height = 300, width = 300)
featplot_gviz("h1SMcG0017844")
dev.off()

## Parenchyma
png("graphics/featplot_scrna_parenchyma.png", height = 300, width = 300)
featplot_gviz("h1SMcG0000267")
dev.off()

## Protonephridia
png("graphics/featplot_scrna_protonephridia.png", height = 300, width = 300)
featplot_gviz("h1SMcG0003681")
dev.off()

## Secretory
png("graphics/featplot_scrna_secretory.png", height = 300, width = 300)
featplot_gviz("h1SMcG0017255")
dev.off()

# Then the scATAC (we run this in the respective markdown/environment with the scatac object loaded)

featplot_atac <- function(x){
  FeaturePlot(
    smed_cisreg_scatac,
    features = x,
    cols = c("#DAE7F2","#174fbc"),
    order = TRUE,
    pt.size = 1.25
    )+NoAxes()+NoLegend()
}

## Early Epid Prog
png("~/projects/smed_cisreg/graphics/featplot_scatac_early_epid_prog.png", height = 300, width = 300)
featplot_atac("chr2-h1-195456153-195457025")
dev.off()

## Late Epid Prog
png("~/projects/smed_cisreg/graphics/featplot_scatac_late_epid_prog.png", height = 300, width = 300)
featplot_atac("chr2-h1-168383475-168384395")
dev.off()

## Epid
# center from beginning of peak to end of gene
png("~/projects/smed_cisreg/graphics/featplot_scatac_epid.png", height = 300, width = 300)
featplot_atac("chr2-h1-261657857-261658707")
dev.off()

## Phag
png("~/projects/smed_cisreg/graphics/featplot_scatac_phag.png", height = 300, width = 300)
featplot_atac("chr2-h1-14407194-14408173")
dev.off()

## Basal/Goblet
png("~/projects/smed_cisreg/graphics/featplot_scatac_basal_goblet.png", height = 300, width = 300)
featplot_atac("chr2-h1-125376387-125377130")
dev.off()

## Musc
png("~/projects/smed_cisreg/graphics/featplot_scatac_musc.png", height = 300, width = 300)
featplot_atac("chr1-h1-113840624-113841542")
dev.off()

## Neur
png("~/projects/smed_cisreg/graphics/featplot_scatac_neur.png", height = 300, width = 300)
featplot_atac("")
FeaturePlot(smed_cisreg_scatac, features = "chr3-h1-57579192-57580102", cols = c("#DAE7F2","#174fbc"), order = TRUE, pt.size = 1.25)+NoAxes()+NoLegend()
dev.off()

## Parenchyma
png("~/projects/smed_cisreg/graphics/featplot_scatac_parenchyma.png", height = 300, width = 300)
featplot_atac("chr1-h1-9570978-9571995")
dev.off()

## Protonephridia
png("~/projects/smed_cisreg/graphics/featplot_scatac_protonephridia.png", height = 300, width = 300)
featplot_atac("chr1-h1-142152671-142153674")
dev.off()

## Secretory
png("~/projects/smed_cisreg/graphics/featplot_scatac_secretory.png", height = 300, width = 300)
featplot_atac("chr3-h1-35797953-35798876")
dev.off()

genes_coords <-
  c(
    "h1SMcG0014249", # EEP
    "h1SMnG0022130", # LEP
    "h1SMcG0016267", # Epd
    "h1SMcG0009317", # Phg
    "h1SMcG0002232", # B_G
    "h1SMcG0003060", # Mus
    "h1SMcG0017844", # Neu
    "h1SMcG0000267", # Par
    "h1SMcG0003681", # PrN
    "h1SMcG0017255"  # Sec
  )


smed_fig1_coords <- read.table("smed_fig1_coords.tsv", sep = "\t", header = TRUE)

pdf("graphics/fig1_coords_scatac_notracks.pdf", height = 4, width = 2)
for(i in 1:nrow(smed_fig1_coords)){
  
  region_coords <- smed_fig1_coords[i,]
  l <- abs(region_coords[1,2]-region_coords[1,3])
  region_coords[1,2] <- region_coords[1,2] - (l*0.5)
  region_coords[1,3] <- region_coords[1,3] + (l*0.5)
  
  # Plotting
  genomeviz(
    coords = region_coords,
    tracks = tracks_table_scatac,
    basic_tracks = NULL,
    y_axis = 100
  )
  
  print(i)
}
dev.off()

pdf("graphics/fig1_coords_scrna_notracks.pdf", height = 4, width = 2)
for(i in 1:nrow(smed_fig1_coords)){
  
  region_coords <- smed_fig1_coords[i,]
  l <- abs(region_coords[1,2]-region_coords[1,3])
  region_coords[1,2] <- region_coords[1,2] - (l*0.5)
  region_coords[1,3] <- region_coords[1,3] + (l*0.5)
  
  # Plotting
  genomeviz(
    coords = region_coords,
    tracks = tracks_table_scrna,
    basic_tracks = NULL,
    y_axis = 300
  )
  
  print(i)
}
dev.off()

pdf("graphics/fig1_coords.pdf", height = 8, width = 2)
for(i in 1:nrow(smed_fig1_coords)){
  
  region_coords <- smed_fig1_coords[i,]
  l <- abs(region_coords[1,2]-region_coords[1,3])
  region_coords[1,2] <- region_coords[1,2] - (l*0.5)
  region_coords[1,3] <- region_coords[1,3] + (l*0.5)
  
  # Plotting
  genomeviz(
    coords = region_coords,
    tracks = tracks_table,
    basic_tracks = list(gm, peaks, peaks_sc),
    y_axis = 100
  )
  
  print(i)
}
dev.off()


