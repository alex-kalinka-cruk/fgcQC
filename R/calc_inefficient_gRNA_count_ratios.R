#' calc_inefficient_gRNA_count_ratios
#'
#' Calculates the ratio of normalized counts for inefficient sgRNA PAM-proximal sequences relative to all guides (taking the median of each).
#'
#' @param counts A data frame of normalized counts for each sample in the study (samples as columns, gRNAs as rows).
#' @param library A data frame containing the library file in which the first column gives the sgRNA sequence and the second column gives the sgRNA ID.
#'
#' @return A data frame containing three columns: `SampleName`, `norm_counts_GCC_ratio` and `norm_counts_TT_ratio`.
#' @author Alex T. Kalinka, \email{alex.kalinka@@cancer.org.uk}
#' @importFrom dplyr mutate rowwise ungroup summarise_if group_by
#' @references Graf, R. et al. (2019) sgRNA sequence motifs blocking efficient CRISPR/Cas9-mediated gene editing. Cell Rep 26: 1098-1103.
#' @export
calc_inefficient_gRNA_count_ratios <- function(counts, library){
  tryCatch({
    library %<>%
      fgcQC::extract_proximal_4bases_PAM()
    counts_median <- apply(counts[,3:ncol(counts)],2,median)
    counts %<>%
      dplyr::mutate(proximal_4bases_pam = library$proximal_4bases_pam[match(sgRNA, library$V2)]) %>%
      dplyr::group_by(proximal_4bases_pam) %>%
      dplyr::summarise_if(is.numeric, median, na.rm = TRUE) %>%
      dplyr::ungroup()
    ret <- data.frame(SampleName = colnames(counts)[2:ncol(counts)], stringsAsFactors = F) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(norm_counts_GCC_ratio = unlist(counts[counts$proximal_4bases_pam == "GCC",colnames(counts) == SampleName])/counts_median[names(counts_median) == SampleName],
                    norm_counts_TT_ratio = unlist(counts[counts$proximal_4bases_pam == "TT",colnames(counts) == SampleName])/counts_median[names(counts_median) == SampleName]) %>%
      dplyr::ungroup()
  },
  error = function(e) stop(paste("unable to calculate inefficient gRNA normalized count ratios:",e))
  )
  return(ret)
}
