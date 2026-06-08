clean_sampletable <- function(s){
  s <- s[s$condition != "",]
  s <-
    s[
      s$ctype %in%
        names(table(s$ctype))[ # those cell types
          table(s$ctype) == # for which we have as many instances
            (length(unique(s$condition))*length(unique(s$replicate))) # as possible combinations of condition x replicates (i.e. no missing conditions/reps)
        ],
    ]
  s$condition <- factor(s$condition, levels = unique(s$condition))
  
  return(s)
}