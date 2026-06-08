#' @param coords: a 1-row data.frame, or list, with chromosome, start, end, and a title. Will be used to get the coordinates of the plotting window
#' @param tracks_table: a dataframe with the dataset in bw, bed, etc format to be displayed. Especially tailored for samples and read count-based stuff.
#' @param basic_tracks: a list of gviz datatracks and objects to be passed and plotted by plotTracks (e.g. gene models, cytoband, etc.)
genomeviz <- function(coords, tracks_table, basic_tracks, y_axis){
  require(Gviz)
  require(GenomicRanges)
  options(ucscChromosomeNames = FALSE)
  
  # Parsing coordinates
  coords <- as.list(coords)
  chromosome <- coords$chrom
  from <- coords$from
  to <- coords$to
  coord_name <- coords$name
  
  # Creating tracks
  message("creating tracks")
  #function for creating tracks based on the different files
  quicktrack <- function(x, type = "l"){
    dt <- 
      DataTrack(
        name = x$name,
        range = x$path,
        type = type,
        col = x$col,
        fill = x$fill,
        background.title = x$col
      )
    return(dt)
  }
  if(!("name" %in% colnames(tracks_table))){
    tracks_table$name <- paste(tracks_table[,1],tracks_table[,2], sep = "_")
  }
  
  #create an empty list and fill it with all the tracks, pass this down to the plotting function below
  tracks <- list()
  for(i in 1:nrow(tracks_table)){
    tracks[[i]] <- quicktrack(tracks_table[i,])
  }
  
  if (!(is.null(basic_tracks))){
    tracklist <- c(basic_tracks, tracks)
  } else {
    tracklist <- tracks
  }
  
  # Plotting
  message("plotting")
  plotTracks(
    tracklist, 
    type = "h",
    window = -1,
    main = coord_name,
    cex.main = 0.5,
    chromosome = chromosome,
    ylim = c(0,y_axis),
    from = from, to = to
  )
}
