# The x argument is a character vector where every element is a human-readable GO term. 
abrevi <- function(x) { 
  abr <- character()
  for (i in 1:length(x)) {
    if(nchar(x[i]) > 40){
      a <- paste(
        abbreviate(
          strsplit(x," ")[[i]],
          minlength = 6
        ),
        collapse=" ") # the result of abbreviate is a vector with all the words. We put them back together using paste(..., collapse=" ")
    } else {
      a <- x[i]
    }
    abr <- c(abr,a) # this adds the newly abbreviated term to the list of terms that will be returned
  }
  return(abr) # the list of abbreviated terms is returned
}