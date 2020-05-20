#' calc_low_zero_count_plasmid_gRNAs
#'
#' Caclculate the percentage of low (<30) and zero count gRNAs in the plasmid sample - these should be rare.
#'
#' @param counts A data frame of counts for each sample in the study (samples as columns, gRNAs as rows).
#' @param plasmid_sample_id A character string naming the plasmid sample ID in the column names of `counts`.
#'
#' @return A data frame containing two colums `percent_low_count_plasmid_gRNAs` and `percent_zero_count_plasmid_gRNAs`.
#' @importFrom dplyr summarise
#' @export
calc_low_zero_count_plasmid_gRNAs <- function(counts, plasmid_sample_id){
  if(!plasmid_sample_id %in% colnames(counts))
    stop(paste("could not find",plasmid_sample_id,"in the counts data frame column names:",colnames(counts)))

  tryCatch({
    plc <- sym(plasmid_sample_id)
    zc <- counts %>%
      dplyr::summarise(total_plasmid_reads = sum(!!plc, na.rm=T),
                       percent_low_count_plasmid_gRNAs = 100*sum(!!plc < 30, na.rm=T),
                       percent_zero_count_plasmid_gRNAs = 100*sum(!!plc == 0, na.rm=T))
  },
  error = function(e) stop(paste("unable to calculate plasmid percent zero count gRNA:",e))
  )
  return(zc)
}
