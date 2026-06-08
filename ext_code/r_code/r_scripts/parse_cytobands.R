library(GenomicRanges)
library(Gviz)
library(vroom)

# following https://academic.oup.com/hmg/article/12/9/1037/629726?login=false, we use AT content as a proxy for cytogenetic band staining.
# make intervals file
#bedtools makewindows -g schMedS3_h1.fa.fai -n 500 -i srcwinnum > schMedS3_h1.intervals.bed
# calculate gc content of each interval
#bedtools nuc -fi schMedS3_h1.fa -bed schMedS3_h1.intervals.bed > gc_content.txt
gcontent <- vroom("/mnt/sda/alberto/projects/smed_cisreg/outputs/genomeviz/gc_content.txt")
plot(density(gcontent$`5_pct_at`))
qs <- quantile(gcontent$`5_pct_at`,probs=c(0,0.5,0.65,0.85,0.95,1))
bands <- gcontent[,c(1:4)]; colnames(bands) <- c("chrom","chromStart","chromEnd","name")
bands$gieStain <- 
  cut(
    gcontent$`5_pct_at`,
    breaks = qs,
    labels = c("gneg","gpos25","gpos50","gpos75","gpos100"),
    include.lowest = TRUE
  )
saveRDS(bands,"~/projects/smed_cisreg/outputs/genomeviz/cytobands.rds")
