#' mock_library_cleanr_tsv
#'
#' Mock a 'cleanr.tsv' gRNA library annotation file, write to a temporary file, and return the temporary file path.
#'
#' @param combined_counts A path to a valid combined counts csv file (produced by the AZ-CRUK CRISPR analysis pipeline).
#'
#' @return A file path to a temporary mock 'cleanr' file.
#' @export
mock_library_cleanr_tsv <- function(combined_counts){
  tryCatch({
    counts <- read.delim(combined_counts, sep="\t", header=T, stringsAsFactors = F)
    library <- data.frame(seq = "AAATTTGGGGCC", CODE = counts[,1], GENES = counts[,2], stringsAsFactors = F)
    if(nrow(library) < 3)
      stop(paste("expecting at least 3 gRNAs in",combined_counts))
    # At least one 'TT' gRNA for downstream QC metrics.
    library$seq[2] <- "AAATTTGGGGTT"
    # At least one 'Other' gRNA seq.
    library$seq[3] <- "AAATTTGGGGAA"
    mock_library_file <- tempfile("library_cleanr.")
    write.table(library, mock_library_file, col.names = T, row.names = F, sep = "\t", quote = F)
  },
  error = function(e) stop(paste("mock_library_cleanr_tsv: unable to mock library 'cleanr' file:",e))
  )
  return(mock_library_file)
}
