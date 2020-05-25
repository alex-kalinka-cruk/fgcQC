#' calc_dcorr_GC_content_counts
#'
#' Calculates the distance correlation between normalized counts and the GC content of sgRNAs in a given library.
#'
#' @param counts A data frame of normalized counts for each sample in the study (samples as columns, gRNAs as rows).
#' @param library A data frame containing the library file in which the first column gives the sgRNA sequence and the second column gives the sgRNA ID.
#'
#' @return A data frame with a `SampleName` column and a `distcorr_GC_content_counts` column.
#' @importFrom dplyr %>% mutate rowwise ungroup filter
#' @importFrom energy dcor
#' @importFrom magrittr %<>%
#' @references Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007), Measuring and Testing Dependence by Correlation of Distances, Annals of Statistics, Vol. 35 No. 6, pp. 2769-2794.
#' @export
calc_dcorr_GC_content_counts <- function(counts, library){
  tryCatch({
    library %<>%
      fgcQC::calc_GC_percent_library() %>%
      dplyr::filter(V2 %in% counts$sgRNA)
    counts %<>%
      dplyr::filter(sgRNA %in% library$V2)
    ret <- data.frame(SampleName = colnames(counts)[3:ncol(counts)]) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(distcorr_GC_content_counts = energy::dcor2d(counts[,colnames(counts) == SampleName],
                                                                library$GC_percent[match(counts$sgRNA, library$V2)],
                                                                type = "U")) %>%
      dplyr::ungroup()
  },
  error = function(e) stop(paste("unable to calculate dist corr for GC vs counts:",e))
  )
  return(ret)
}
