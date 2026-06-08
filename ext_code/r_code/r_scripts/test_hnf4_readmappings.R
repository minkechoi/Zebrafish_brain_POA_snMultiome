## Genome Visualisation
library(Gviz)
library(GenomicRanges)
library(rtracklayer)
source("~/projects/smed_cisreg/code/r_code/functions/genomeviz.R")
options(ucscChromosomeNames=FALSE)

# Load the genomic coordinates
coords_of_interest <-
  data.frame(
    chrom = "chr3_h1",
    from = 107846166,
    to = 107855701,
    name = "h1SMcG0019688 / Hnf4"
  )

# Load the track information table
tracks_table <- read.delim2(
  file = "~/projects/smed_cisreg/outputs/hnf_KD/bams_by_experiment/tracks.tsv",
  sep = "\t",
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


dsrna_bed <- import.bed("/mnt/sda/alberto/projects/smed_cisreg/outputs/hnf_KD/dsrna.bed")
hnf4_dsrna <- AnnotationTrack(range=dsrna_bed, genome="Smed", name = "dsRNA", fill = "red", background.title = "darkgrey")

ccc <- coords_of_interest
l <- abs(ccc[1,2]-ccc[1,3])
ccc[1,2] <- ccc[1,2] - (l*0.5)
ccc[1,3] <- ccc[1,3] + (l*0.5)

# Plotting
genomeviz(
  coords = ccc,
  tracks = tracks_table[3:6,],
  basic_tracks = list(itrack, axisTrack, hnf4_dsrna,gm),
  y_axis = 20
)

