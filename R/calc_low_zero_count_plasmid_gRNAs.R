# Helper function.
.median <- function(...){
  x <- c(unlist(list(...)))
  return(median(x))
}


#' calc_low_zero_count_plasmid_gRNAs
#'
#' Caclculate the percentage of low (<30) and zero count gRNAs in the plasmid sample - these should be rare.
#'
#' @param counts A data frame of counts for each sample in the study (samples as columns, gRNAs as rows).
#' @param plasmid_sample_id A character string naming the plasmid sample ID in the column names of `counts`. If several plasmid samples are given, the median count will be used.
#'
#' @return A data frame containing two colums `percent_low_count_plasmid_gRNAs` and `percent_zero_count_plasmid_gRNAs`.
#' @importFrom dplyr summarise mutate
#' @importFrom purrr pmap_dbl
#' @export
calc_low_zero_count_plasmid_gRNAs <- function(counts, plasmid_sample_id){
  if(length(setdiff(plasmid_sample_id, colnames(counts))) > 0)
    stop(paste("could not find",plasmid_sample_id,"in the counts data frame column names:",colnames(counts)))

  tryCatch({
    zc <- as.data.frame(counts[,plasmid_sample_id]) %>%
      dplyr::mutate(median_count = purrr::pmap_dbl(
        .l = .,
        .f = .median
      )) %>%
      dplyr::summarise(percent_low_count_plasmid_gRNAs = 100*sum(median_count < 30, na.rm=T)/n(),
                       percent_zero_count_plasmid_gRNAs = 100*sum(median_count == 0, na.rm=T)/n())
  },
  error = function(e) stop(paste("unable to calculate plasmid percent zero count gRNA:",e))
  )
  return(zc)
}
