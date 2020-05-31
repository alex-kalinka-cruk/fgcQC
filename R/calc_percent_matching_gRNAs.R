#' calc_percent_matching_gRNAs
#'
#' Calculates the percentage of reads in each sample that can be assigned to any of the gRNA sequences.
#'
#' @param b2f_samples A data frame of sample information derived from `bcl2fastq2` output - see `fgcQC::extract_b2f_json`. The `SampleName` column should have been replaced with sample names hailing from the analysis config.
#' @param counts A data frame of counts for each sample in the study (samples as columns, gRNAs as rows).
#'
#' @return A data frame containing per-sample percent matching gRNA values in a column called `percent_reads_matching_gRNAs`.
#' @author Alex T. Kalinka, \email{alex.kalinka@@cancer.org.uk}
#'
#' @importFrom dplyr mutate
#' @export
calc_percent_matching_gRNAs <- function(b2f_samples, counts){
  tryCatch({
    totals <- colSums(counts[,3:ncol(counts)])
    pm <- b2f_samples %>%
      dplyr::mutate(percent_reads_matching_gRNAs = 100*totals[match(SampleName, names(totals))]/NumberReads)
  },
  error = function(e) paste(stop("calc_percent_matching_gRNAs: unable to calculate reads matching gRNA percent:",e))
  )
  return(pm)
}
