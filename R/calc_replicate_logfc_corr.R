#' calc_replicate_logfc_corr
#'
#' Calculates logFC replicate Pearson correlations in accordance with recommendations by Hanna & Doench (2020).
#'
#' @param logfc A data frame with columns `log2FC.repl_1` and `log2FC.repl_2`.
#' @param suffix A character string naming the comparison: one of `ctrl_plasmid` and `treat_plasmid`.
#'
#' @importFrom rlang :=
#' @importFrom tibble tibble
#' @references Hanna, R. E. and Doench, J. G. 2020. Design and analysis of CRISPR-Cas experiments. Nature Biotechnology.
#' @export
calc_replicate_logfc_corr <- function(logfc, suffix){
  tryCatch({
    if(!is.null(logfc)){
      ret <- tibble::tibble(!!paste("repl_log2FC_pearson_corr", suffix, sep=".") := cor(logfc$log2FC.repl_1, logfc$log2FC.repl_2,
                                                                                      use = "pairwise.complete.obs"))
    }else{
      ret <- tibble::tibble(!!paste("repl_log2FC_pearson_corr", suffix, sep=".") := NA)
    }
    return(ret)
  },
  error = function(e) stop(paste("calc_replicate_logfc_corr: unable to calculate replicate logFC Pearson correlation:",e))
  )
}
