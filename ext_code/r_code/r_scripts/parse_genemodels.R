library(GenomicRanges)
library(seqQTL)
rosetta <- read.delim2("/mnt/sda/Standard_References/Smed_Rink/20231127_Rosetta.tsv")


#' Load gene models in proper format
#' We need to provide the info on strand, exons from same manuscript, 
#' and likely show just the gene models and not the transcript models.
#' For this we will load the GTF as a GRanges object that can be properly parsed.
#' The idea is: making a GRanges object with the same features as the geneModels object from the tutorial of GViz
smed_gtf <- rtracklayer::import("~/projects/smed_cisreg/data/standard_references/schMedS3_h1.gtf")

# subset only exons and utrs
my_gr <- split(smed_gtf,smed_gtf$type)
my_gr <- c(my_gr[["exon"]],my_gr[["five_prime_UTR"]],my_gr[["three_prime_UTR"]])

# for every exon and utr, fix the transcript and gene using the gene_id and transcript_id
mcols(my_gr, level="within")[, "feature"] <- mcols(my_gr, level="within")[, "type"]
mcols(my_gr, level="within")[, "gene"] <- mcols(my_gr, level="within")[, "gene_id"]
mcols(my_gr, level="within")[, "transcript"] <- mcols(my_gr, level="within")[, "transcript_id"]
mcols(my_gr, level="within")[, "exon"] <- mcols(my_gr, level="within")[, "ID"]
mcols(my_gr, level="within")[, "symbol"] <- mcols(my_gr, level="within")[, "gene"]

my_gr$feature <- droplevels(my_gr$feature)
levels(my_gr$feature) <- c("protein-coding","utr5","utr3")

# Fix the utrs
# first we split the granges by type of feature
my_gr <- split(my_gr,my_gr$feature)

# in utr5, exon column must be the name of the closest exon of the same gene from that utr5, that is DOWNstream of the utr5 
#Find the nearest downstream exon range for each utr5 range
nearest_exons <- follow(my_gr$utr5, my_gr$`protein-coding`)

# Extract the "ID" of the nearest downstream exon range
nearest_exon_ids <- my_gr$`protein-coding`$ID[nearest_exons]

# Add the "nearest_exon_id" column to the utr5 ranges
my_gr$utr5$exon <- nearest_exon_ids

# in utr3, exon column must be the name of the closest exon of the same gene from that utr3, that is UPstream of the utr3 
#Find the nearest downstream exon range for each utr5 range
nearest_exons <- precede(my_gr$utr3, my_gr$`protein-coding`)

# Extract the "ID" of the nearest downstream exon range
nearest_exon_ids <- my_gr$`protein-coding`$ID[nearest_exons]

# Add the "nearest_exon_id" column to the utr5 ranges
my_gr$utr3$exon <- nearest_exon_ids

# Add an "exon" column to, well, exons, since column "ID" will disappear
my_gr$`protein-coding`$exon <- my_gr$`protein-coding`$ID

# Reassemble the granges after working out the exon IDs
my_gr <- c(my_gr[["protein-coding"]],my_gr[["utr5"]],my_gr[["utr3"]])

# keep good columns only
genemodels <- my_gr[,c("feature","gene","exon","transcript", "symbol")] #can we save this as an object we can load? maybe a .rds?

# Save as object to load in the markdown
saveRDS(genemodels,"~/projects/smed_cisreg/outputs/genomeviz/smed_cisreg_genomeviz_genemodels.rds")

# Create an object with longest isoforms only
smed_genemodels_longestIso <- genemodels[genemodels$transcript %in% rosetta$longest_isoform]

# Save the object with longest isoforms only as a GFF3
GR2gff(smed_genemodels_longestIso, "~/projects/smed_cisreg/outputs/genomeviz/smed_cisreg_gff_longestIso.gff3", src = "GenomicRanges", score = ".", phase = ".")

