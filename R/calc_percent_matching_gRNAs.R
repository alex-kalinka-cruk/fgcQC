#' calc_percent_matching_gRNAs
#'
#' Calculates the percentage of reads in each sample that can be assigned to any of the gRNA sequences.
#'
#' @param b2f_samples A data frame of sample information derived from `bcl2fastq2` output - see `extract_b2f_json`. The `SampleName` column should have been replaced with sample names hailing from the analysis config.
#' @param counts A data frame of counts for each sample in the study (samples as columns, gRNAs as rows).
#'
#' @return A data frame containing per-sample percent matching gRNA values in a column called `percent_reads_matching_gRNAs`.
#'
#' @importFrom dplyr mutate
#' @export
calc_percent_matching_gRNAs <- function(b2f_samples, counts){
  if(any(b2f_samples$SampleId == b2f_samples$SampleName))
    stop("the 'SampleName' column in 'b2f_samples' should contain sample names hailing from the analysis config")

  tryCatch({
    totals <- colSums(counts[,3:ncol(counts)])
    pm <- b2f_samples %>%
      dplyr::mutate(percent_reads_matching_gRNAs = 100*NumberReads/totals[match(SampleName, names(totals))])
  },
  error = function(e) paste(stop("unable to calculate matching gRNA percent:",e))
  )
  return(pm)
}
