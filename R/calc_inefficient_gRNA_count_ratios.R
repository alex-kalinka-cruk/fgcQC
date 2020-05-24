#' calc_inefficient_gRNA_count_ratios
#'
#' Calculates the ratio of normalized counts for inefficient sgRNA PAM-proximal sequences relative to all guides (taking the median of each).
#'
#' @param counts A data frame of normalized counts for each sample in the study (samples as columns, gRNAs as rows).
#' @param library A data frame containing the library file in which the first column gives the sgRNA sequence and the second column gives the sgRNA ID.
#'
#' @return A data frame containing three columns: `SampleName`, `norm_counts_GCC_ratio` and `norm_counts_TT_ratio`.
#' @importFrom dplyr mutate filter rowwise
#' @references Graf, R. et al. (2019) sgRNA sequence motifs blocking efficient CRISPR/Cas9-mediated gene editing. Cell Rep 26: 1098-1103.
#' @export
calc_inefficient_gRNA_count_ratios <- function(counts, library){
  tryCatch({
    library %<>%
      extract_proximal_4bases_PAM() %>%
      dplyr::filter(V2 %in% counts$sgRNA)
    counts %<>%
      dplyr::filter(sgRNA %in% library$V2) %>%
      dplyr::mutate(proximal_4bases_pam = library$proximal_4bases_pam[match(sgRNA, V2)])
    ret <- data.frame(SampleName = colnames(counts)[3:ncol(counts)]) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(norm_counts_GCC_ratio = median(counts[counts$proximal_4bases_pam == "GCC",colnames(counts) == SampleName],na.rm=T)/median(counts[,colnames(counts) == SampleName],na.rm=T),
                    norm_counts_TT_ratio = median(counts[counts$proximal_4bases_pam == "TT",colnames(counts) == SampleName],na.rm=T)/median(counts[,colnames(counts) == SampleName],na.rm=T))
  },
  error = function(e) stop(paste("unable to calculate inefficient gRNA count ratios:",e))
  )
  return(ret)
}
