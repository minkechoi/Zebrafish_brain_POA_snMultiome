# Load DB of all jaspar motifs used by ANANSE
ananse_jaspar2020_db <- read.delim2(
  "/mnt/sda/alberto/programs/miniconda3/envs/ananse_venv/lib/python3.10/site-packages/data/motif_databases/JASPAR2020.motif2factors.txt",
  sep = "\t",
  header = TRUE
)
# ANANSE has the motif and the TF ids joined in a single column separated by "_"
# Add a column of just the motif ID
ananse_jaspar2020_db$MAid <- gsub("_.*","",ananse_jaspar2020_db$Motif)

# Load motif2factors output
smed_m2f <-read.delim2("~/projects/smed_cisreg/outputs/ananse/m2f/smed_motif2factors_more_refs_JASPAR_2020/Smed.JASPAR2020.motif2factors.txt")

# some gene IDs are carrying over the ">" symbol, for some reason
smed_m2f$Factor <- sub(">","",smed_m2f$Factor)

# Load profile inference results
smed_proInf_results <- read.delim2("~/projects/smed_cisreg/outputs/jaspar_inference/inferred_matrix.tsv")

# we remove dupped ids separated by "|"
smed_proInf_results$Query <- gsub("\\|.*","",smed_proInf_results$Query)

# we merge with the jaspar db so we can give these motif predictions a matrix_TF id which will facilitate the joining of the tables
smed_proInf_results_m2f <-
  merge(
    smed_proInf_results,
    ananse_jaspar2020_db,
    by.x = "TF.Matrix",
    by.y = "MAid",
    all.x = TRUE
  )

# profile inference was done using JASPAR v2024 matrices and m2f is from 2020. Some matrices are not in the m2f file. We miss 5 TFs this way. But it should be ok because we are also doing gimme's m2f.
# for these tfs, "Y" means "taken from profile inference tool"
smed_inference <- smed_proInf_results_m2f[complete.cases(smed_proInf_results_m2f),][,c(8,2,10,11)] 

# Give it the same column names as m2f output to allow joining the tables
colnames(smed_inference) <- colnames(smed_m2f)

smed_m2f_merged <-
  rbind(
    smed_m2f,
    smed_inference
  )

# which of these motifs have now found an "ortholog"? (predicted by JASPAR's inference). We'll have to remove these entries
filt_out <- 
  which(
    smed_m2f_merged$Motif %in% smed_m2f_merged$Motif[smed_m2f_merged$Factor != "NO ORTHOLOGS FOUND"] & # which row has a motif with an assigned ortholog
      smed_m2f_merged$Factor == "NO ORTHOLOGS FOUND" # and yet it says in this row "NO ORTHOLOGS FOUND"
  )

smed_m2f_merged <- smed_m2f_merged[-filt_out,]

# Load the tfs with an assigned JASPAR matrix by Neiro et al., 2022
smed_tfsneiro <- tfs_all[!(tfs_all$Jaspar.MatrixID %in% c(".","-", "0")),c(1,21)]

# we merge with the jaspar db so we can, again, give them a matrix_TF id which will facilitate the joining of the tables
smed_tfsneiro_m2f <- 
  merge(
    smed_tfsneiro,
    ananse_jaspar2020_db,
    by.x = 2,
    by.y = 5,
    all.x = TRUE
  )

smed_tfsneiro_m2f <- smed_tfsneiro_m2f[complete.cases(smed_tfsneiro_m2f),] # 201 out of 286... these motif matrices also come from JASPAR 2020 supposedly. Where are the missing matrices?

# again we rearrange and tidy to do the joining of tables
smed_tfsneiro_m2f <- smed_tfsneiro_m2f[,c(3,2,5,6)]
colnames(smed_tfsneiro_m2f) <- colnames(ananse_jaspar2020_db)[1:4]

smed_m2f_merged <-
  rbind(
    smed_m2f_merged,
    smed_tfsneiro
  )
# which of these motifs have now found an "ortholog"? (predicted by Neiro et al. when they used JASPAR profile inference)
filt_out <- 
  which(
    smed_m2f_merged$Motif %in% smed_m2f_merged$Motif[smed_m2f_merged$Factor != "NO ORTHOLOGS FOUND"] & # which row has a motif with an assigned ortholog
      smed_m2f_merged$Factor == "NO ORTHOLOGS FOUND" # and yet it says in this row "NO ORTHOLOGS FOUND"
  )

smed_m2f_merged <- smed_m2f_merged[-filt_out,]

# load motif2factors from a previous run to do some comparisons
m2f_previous <- read.delim2("~/projects/smed_cisreg/outputs/ananse/m2f/_obs/smed_motif2factors/Smed.gimme.vertebrate.v5.0.motif2factors.txt")
length(which(unique(m2f_previous$Factor) %in% tfs_all$gene)) / nrow(tfs_all)

length(which(unique(smed_m2f_merged$Factor) %in% tfs_all$gene)) # 401 out of 665. Previously we had 196 out of 665
length(which(unique(smed_m2f_merged$Factor) %in% tfs_all$gene)) / nrow(tfs_all)

# should I just keep the ones we got in the curated TF list or should I just add all? I think I should keep the ones we have in our curated list. Should I have used the same roster of species in the TFanimalDB evidence as here with ANANSE?

write.table(
  smed_m2f_merged,
  "outputs/ananse/m2f/smed_motif2factors_more_refs_JASPAR_2020/merge_m2f_inference_neiroetal/Smed.JASPAR2020.motif2factors.txt",
  sep = "\t",
  quote = F,
  row.names = F
)

smed_m2f_merged_401 <-
  smed_tfs_m2f

smed_m2f_merged_401$Factor <- 
  sapply(
    smed_m2f_merged_401$Factor,
    function(x){
      if(x %in% tfs_all$gene){
        y = x
      } else {
        y = "NO ORTHOLOGS FOUND"
      }
      return(y)
      }
  )

smed_m2f_merged_401 <- unique(smed_m2f_merged_401)

write.table(
  smed_m2f_merged_401,
  "outputs/ananse/m2f/smed_motif2factors_more_refs_JASPAR_2020/merge_m2f_inference_neiroetal/core_common_tfs/Smed.JASPAR2020.motif2factors.txt",
  sep = "\t",
  quote = F,
  row.names = F
)


