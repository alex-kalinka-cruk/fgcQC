#' calc_log2_fold_change_gRNAs
#'
#' Calculates log2 fold changes for sgRNAs from normalized counts. Median counts are taken across replicates.
#'
#' @param counts A data frame of normalized sgRNA counts for each sample in the study (samples as columns, gRNAs as rows).
#' @param ref A character vector naming the reference sample(s).
#' @param comp A character vector naming the comparison sample(s) to be contrasted against `ref` samples.
#' @param pseudo_count A small integer specifying a pseudo-count to add to the counts to avoid log2(0) cases. Defaults to 5.
#'
#' @return A data frame of log2 fold changes in which the columns are samples and the rows are sgRNAs
#'
#' @importFrom dplyr rowwise mutate ungroup group_by
#' @importFrom magrittr %<>%
#' @export
calc_log2_fold_change_gRNAs <- function(counts, ref, comp, pseudo_count = 5){
  if(!all(ref %in% colnames(counts)))
    stop(paste("unable to find",ref,"in column names of 'counts'"))
  if(!all(comp %in% colnames(counts)))
    stop(paste("unable to find",comp,"in column names of 'counts'"))

  tryCatch({
   counts$ref_median <- apply(counts[,ref], 1, median, na.rm = T)
   counts$comp_median <- apply(counts[,comp], 1, median, na.rm = T)
   counts %<>%
     dplyr::rowwise() %>%
     dplyr::mutate(log2FC = log2(ref_median + pseudo_count) - log2(comp_median + pseudo_count)) %>%
     dplyr::ungroup()
  },
  error = function(e) stop(paste("unable to calculate log2 fold changes at the gene level:",e))
  )
  return(counts)
}
