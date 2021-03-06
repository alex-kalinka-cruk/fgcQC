#' calc_inefficient_gRNA_logfc
#'
#' Calculates the logFC difference for inefficient sgRNA PAM-proximal sequences relative to all guides (taking the median of each).
#'
#' @param logfc A data frame of log2 fold change data for each sample in the study (samples as columns, gRNAs as rows) as generated by `fgcQC::calc_log2_fold_change_gRNAs`.
#' @param library A data frame containing the library file in which the first column gives the sgRNA sequence and the second column gives the sgRNA ID.
#' @param col_suffix A character string providing a name to go at the end of the output column name (e.g. `ctrl_plasmid`).
#'
#' @return A data frame containing two columns: `log2FC_GCC_diff` and `log2FC_TT_diff`.
#' @author Alex T. Kalinka, \email{alex.kalinka@@cancer.org.uk}
#' @importFrom dplyr mutate filter rowwise sym
#' @importFrom rlang :=
#' @importFrom tibble tibble
#' @references Graf, R. et al. (2019) sgRNA sequence motifs blocking efficient CRISPR/Cas9-mediated gene editing. Cell Rep 26: 1098-1103.
#' @export
calc_inefficient_gRNA_logfc <- function(logfc, library, col_suffix){
  tryCatch({
    gcc_col <- dplyr::sym(paste("log2FC_GCC_diff", col_suffix, sep="."))
    tt_col <- dplyr::sym(paste("log2FC_TT_diff", col_suffix, sep="."))
    library %<>%
      fgcQC::extract_proximal_4bases_PAM() %>%
      dplyr::filter(V2 %in% logfc$sgRNA)
    logfc %<>%
      dplyr::filter(sgRNA %in% library$V2) %>%
      dplyr::mutate(proximal_4bases_pam = library$proximal_4bases_pam[match(sgRNA, library$V2)])
    ret <- tibble::tibble(!!gcc_col := median(logfc$log2FC[logfc$proximal_4bases_pam == "GCC"], na.rm=T) - median(logfc$log2FC, na.rm=T),
                          !!tt_col := median(logfc$log2FC[logfc$proximal_4bases_pam == "TT"], na.rm=T) - median(logfc$log2FC, na.rm=T))
  },
  error = function(e) stop(paste("unable to calculate inefficient gRNA logFC ratios:",e))
  )
  return(ret)
}
